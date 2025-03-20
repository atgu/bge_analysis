from typing import List, Optional

import hailtop.batch as hb
from hailtop.batch.job import Job
import hailtop.fs as hfs

from ..globals import Chunk, SampleGroup, get_bucket


def copy_temp_crams_job(b: hb.Batch,
                        jg: hb.JobGroup,
                        sample_group: SampleGroup,
                        start: int,
                        end: int,
                        cpu: int,
                        memory: str) -> Job:
    # this job is idempotent because we use the -n or no-clobber option

    copy_list = sample_group.write_gcloud_cram_copy_list(start, end)
    copy_list_input = b.read_input(copy_list)

    j = jg.new_bash_job(name=f'copy/sample-group-{sample_group.sample_group_index}/{start}-{end}',
                        attributes={'task': 'copy'})
    j.image('google/cloud-sdk:512.0.0-slim')

    copy = f'cat {copy_list_input} | gcloud storage cp -n -I {sample_group.remote_cram_temp_dir}'
    j.command(f'{copy} || (sleep 2 && {copy}) || (sleep 4 && {copy}) || (sleep 8 && {copy}) || (sleep 16 && {copy})')
    j.cpu(cpu)
    j.memory(memory)

    return j


def phase(b: hb.Batch,
          jg: hb.JobGroup,
          output_file: str,
          phase_file_exists: bool,
          sample_group: SampleGroup,
          chunk: Chunk,
          cram_remote_tmp_path: str,
          mount_point: str,
          crams_list: hb.ResourceFile,
          glimpse_remote_checkpoint_file: str,
          fasta: hb.ResourceGroup,
          sample_ploidy_list: hb.ResourceFile,
          docker: str,
          cpu: int,
          memory: str,
          use_checkpoint: bool,
          impute_reference_only_variants: bool,
          call_indels: bool,
          n_burn_in: Optional[int],
          n_main: Optional[int],
          effective_population_size: Optional[int]) -> Optional[Job]:
    sample_group_index = sample_group.sample_group_index

    glimpse_checkpoint_file_input = b.read_input(glimpse_remote_checkpoint_file)

    if use_checkpoint and phase_file_exists:
        return None

    j = jg.new_bash_job(name=f'phase/sample-group-{sample_group_index}/{chunk.chunk_contig}/{chunk.chunk_idx}',
                        attributes={'sample-group-index': str(sample_group_index),
                                    'contig': str(chunk.chunk_contig),
                                    'chunk-index': str(chunk.chunk_idx),
                                    'task': 'phase'})

    j.image(docker)
    j.storage('20Gi')
    j.cpu(cpu)
    j.memory(memory)

    cram_bucket = get_bucket(cram_remote_tmp_path)
    j.cloudfuse(cram_bucket, mount_point, read_only=True)

    chunk_input = b.read_input(chunk.path)

    j.declare_resource_group(phased={'bcf': '{root}.bcf',
                                     'csi': '{root}.bcf.csi',
                                     'coverage_metrics': '{root}_stats_coverage.txt.gz'})

    extra_args = ''
    if impute_reference_only_variants:
        extra_args += ' --impute-reference-only-variants'
    if call_indels:
        extra_args += ' --call-indels'
    if n_burn_in is not None:
        extra_args += f' --burnin {n_burn_in}'
    if n_main is not None:
        extra_args += f' --main {n_main}'
    if effective_population_size is not None:
        extra_args += f' --ne {effective_population_size}'

    if chunk.is_non_par:
        extra_args += f' --samples-file {sample_ploidy_list}'

    phase_cmd = f'''
set -e

while true; do
  gcloud storage cp checkpoint.bin {glimpse_remote_checkpoint_file} || true;
  sleep 60;
done &

cmd="/bin/GLIMPSE2_phase \
    --reference {chunk_input} \
    --output {j.phased.bcf} \
    --threads {cpu} \
    --bam-list {crams_list} \
    --fasta {fasta.fasta} \
    --checkpoint-file-out checkpoint.bin {extra_args}"

if [ -s {glimpse_checkpoint_file_input} ]; then
    cmd="$cmd --checkpoint-file-in {glimpse_checkpoint_file_input}" 
fi

#check for read error which corresponds exactly to end of cram/bam block.  
#This currently triggers a warning message from htslib, but doesn't return any error.
#We need to make sure that stderr is maintained since cromwell looks for oom strings
#in stderr

eval $cmd 2> >(tee glimpse_stderr.log >&2) 

if grep -q "EOF marker is absent" glimpse_stderr.log; then 
    echo "An input file appears to be truncated.  This may be either a truly truncated file which needs to be fixed, or a networking error which can just be retried."
    exit 1
fi

touch {j.phased.csi}
touch {j.phased.coverage_metrics}
'''

    j.command(phase_cmd)

    b.write_output(j.phased.bcf, output_file + '.bcf')
    b.write_output(j.phased.csi, output_file + '.bcf.csi')
    b.write_output(j.phased.coverage_metrics, output_file + '_stats_coverage.txt.gz')

    return j


def ligate(b: hb.Batch,
           jg: hb.JobGroup,
           sample_group: SampleGroup,
           contig: str,
           docker: str,
           cpu: int,
           memory: str,
           storage: str,
           chunk_outputs: List[hb.ResourceGroup],
           output_file: str,
           ref_dict: hb.ResourceFile,
           use_checkpoint: bool) -> Optional[Job]:
    if use_checkpoint and hfs.exists(output_file + '.vcf.bgz'):
        return None

    j = jg.new_bash_job(name=f'ligate/sample-group-{sample_group.sample_group_index}/{contig}',
                        attributes={'task': 'ligate'})
    j.image(docker)
    j.cpu(cpu)
    j.memory(memory)
    j.storage(storage)

    j.declare_resource_group(ligated={'vcf': '{root}.vcf.gz',
                                      'tbi': '{root}.vcf.gz.tbi',
                                      'md5sum': '{root}.md5sum'})

    chunk_files_str = '\n'.join([str(chunk.bcf) for chunk in chunk_outputs])

    touch_files_str = '\n'.join([f'touch {chunk.csi}' for chunk in chunk_outputs])

    ligate_cmd = f'''
cat > input_list.txt <<EOF
{chunk_files_str}
EOF

cat > touch.sh <<EOF
{touch_files_str}
EOF

sh touch.sh

/bin/GLIMPSE2_ligate --input input_list.txt --output ligated.vcf.gz --threads {cpu}

# Set correct reference dictionary
bcftools view -h --no-version ligated.vcf.gz > old_header.vcf        
java -jar /picard.jar UpdateVcfSequenceDictionary -I old_header.vcf --SD {ref_dict} -O new_header.vcf        
bcftools reheader -h new_header.vcf -o {j.ligated.vcf} ligated.vcf.gz

tabix {j.ligated.vcf}
md5sum {j.ligated.vcf} | awk '{{ print $1 }}' > {j.ligated.md5sum}
touch {j.ligated.tbi}  # this is a dummy operation; todo to figure out why this is necessary
'''

    j.command(ligate_cmd)

    b.write_output(j.ligated.vcf, output_file + '.vcf.bgz')
    b.write_output(j.ligated.tbi, output_file + '.vcf.bgz.tbi')
    b.write_output(j.ligated.md5sum, output_file.replace('.vcf.bgz', '') + '.md5sum')

    return j


def delete_temp_files_job(b: hb.Batch,
                          jg: hb.JobGroup,
                          sample_group: SampleGroup,
                          save_checkpoints: bool) -> Job:
    j = jg.new_bash_job(attributes={'name': f'delete/sample-group-{sample_group.sample_group_index}',
                                    'task': 'delete'})
    j.cpu(1)
    j.image('google/cloud-sdk:512.0.0-slim')

    j.command('set +e')

    # always delete temp cram files as these are expensive to store!
    j.command(f'gcloud storage rm {sample_group.remote_cram_temp_dir}*')

    if not save_checkpoints:
        j.command(f'gcloud storage rm {sample_group.temp_dir}/copy-crams/*/*')
        # FIXME: these need to change if we're saving the coverage_stats
        j.command(f'gcloud storage rm {sample_group.temp_dir}/phase/*/*/*')
        j.command(f'gcloud storage rm {sample_group.temp_dir}/phase/*/*')
        j.command(f'gcloud storage rm {sample_group.temp_dir}/phase/*')
        j.command(f'gcloud storage rm {sample_group.temp_dir}/ligate/*')

    return j


def write_success(b: hb.Batch, jg: hb.JobGroup, sample_group: SampleGroup, docker: str) -> Job:
    j = jg.new_bash_job(attributes={'name': f'write-success/sample-group-{sample_group.sample_group_index}',
                                    'task': 'write-success'})
    j.cpu(1)
    j.image(docker)
    j.command(f'''
python3 << EOF
import hailtop.fs as hfs
with hfs.open("{sample_group.success_file}", "w") as f:
    f.write("0")
EOF
''')

    return j


def union_sample_groups_from_vcfs(b: hb.Batch,
                                  jg: hb.JobGroup,
                                  vcf_paths: hb.ResourceFile,
                                  output_path: str,
                                  docker: str,
                                  cpu: int,
                                  memory: str,
                                  storage: str,
                                  billing_project: str,
                                  remote_tmpdir: str,
                                  regions: str,
                                  use_checkpoints: bool,
                                  contig: str) -> Optional[Job]:
    if use_checkpoints:
        if output_path.endswith('.vcf.bgz') and hfs.exists(output_path):
            return None
        if output_path.endswith('.mt') and hfs.exists(output_path + '/_SUCCESS'):
            return None

    j = jg.new_bash_job(attributes={'name': f'union/{contig}'})
    j.cpu(cpu)
    j.image(docker)
    j.storage(storage)
    j.memory(memory)
    j.spot(False)


    def add_info_if_needed(mt):
        return f'{mt} = {mt}.annotate_rows(info={mt}.info.annotate(INFO={mt}.info.get("INFO", hl.null(hl.tarray(hl.tfloat64)))))'


    cmd = f"""
hailctl config set batch/billing_project "{billing_project}"
hailctl config set batch/remote_tmpdir "{remote_tmpdir}"
hailctl config set batch/regions "{','.join(regions)}"

python3 << EOF
import hail as hl
import hailtop.fs as hfs
import os
from typing import List
import pandas as pd

batch_id = int(os.environ['HAIL_BATCH_ID'])

hl.init(backend="batch", app_name="union-{contig}")  # FIXME: batch_id=batch_id

paths = []
sample_sizes = []
with open("{vcf_paths}", 'r') as f:
    for line in f:
        path, sample_size = line.rstrip("\\n").split('\\t')
        paths.append(path)
        sample_sizes.append(int(sample_size))

mt_left = hl.import_vcf(paths[0], reference_genome="GRCh38")
{add_info_if_needed('mt_left')}
mt_left = mt_left.annotate_rows(info=mt_left.info.annotate(N=sample_sizes[0], AF=mt_left.info.AF[0], INFO=mt_left.info.INFO[0], RAF=mt_left.info.RAF[0]))
mt_left = mt_left.annotate_rows(**{{"info_0": mt_left.info}})

for idx, path in enumerate(paths[1:]):
    mt_right = hl.import_vcf(path, reference_genome="GRCh38")
    {add_info_if_needed('mt_right')}
    mt_right = mt_right.annotate_rows(info=mt_right.info.annotate(N=sample_sizes[idx], AF=mt_right.info.AF[0], INFO=mt_right.info.INFO[0], RAF=mt_right.info.RAF[0]))
    mt_left = mt_left.union_cols(mt_right,
                                 drop_right_row_fields=False,
                                 row_join_type='outer')

mt = mt_left

n_samples = mt.count_cols()
n_batches = len(paths)

mt = mt.annotate_rows(info=mt.info.annotate(AF=hl.array([mt[f"info_{{i}}"].AF for i in range(n_batches)])))
mt = mt.annotate_rows(info=mt.info.annotate(INFO=hl.array([mt[f"info_{{i}}"].INFO for i in range(n_batches)])))
mt = mt.annotate_rows(info=mt.info.annotate(N=hl.array([mt[f"info_{{i}}"].N for i in range(n_batches)])))


def GLIMPSE_AF(mt):
    return hl.sum(hl.map(lambda af, n: af * n, mt.info.AF, mt.info.N)) / n_samples


def GLIMPSE_INFO(mt):
    return hl.if_else((GLIMPSE_AF(mt) == 0) | (GLIMPSE_AF(mt) == 1),
                      1,
                      1 - hl.sum(hl.map(lambda af, n, info: (1 - info) * 2 * n * af * (1 - af), mt.info.AF, mt.info.N, mt.info.INFO)) / (2 * n_samples * GLIMPSE_AF(mt) * (1 - GLIMPSE_AF(mt))))


mt = mt.annotate_rows(info=mt.info.annotate(AF=GLIMPSE_AF(mt), INFO=GLIMPSE_INFO(mt)))

mt = mt.annotate_rows(info=mt.info.drop('N'))
mt = mt.drop(*[f'info_{{i}}' for i in range(n_batches)])
mt = mt.drop(*[f'rsid_{{i}}' for i in range(1, n_batches)])
mt = mt.drop(*[f'qual_{{i}}' for i in range(1, n_batches)])
mt = mt.drop(*[f'filters_{{i}}' for i in range(1, n_batches)])

output_path = "{output_path}"

if output_path.endswith('.mt'):
    mt.write(output_path)
    mt_count = hl.read_matrix_table(output_path)
    print(mt_count.count())
    print(mt.describe())
else:
    assert output_path.endswith('.vcf.bgz')
    hl.export_vcf(mt, output_path, tabix=True)
EOF
"""

    j.command(cmd)

    return j
