import copy
from typing import List, Optional

import hailtop.batch as hb
import hailtop.batch_client.aioclient as bc
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


possible_resources = []
for cpu in (1, 2, 4, 8, 16):
    for memory in ('lowmem', 'standard', 'highmem'):
        possible_resources.append((cpu, memory))


class JobInfo:
    @staticmethod
    def from_json(d: dict):
        return JobInfo(d['name'], d['batch_id'], d['job_id'], d['state'], exit_code=d['exit_code'])

    def __init__(self, name, batch_id, job_id, state, exit_code):
        self.name = name
        self.batch_id = batch_id
        self.job_id = job_id
        self.state = state
        self.exit_code = exit_code
        self.attempt_number = None

    def is_oom(self):
        return self.exit_code is not None and self.exit_code == 137

    def should_be_rerun(self):
        return self.is_oom() or self.state in ('Failed', 'Error')

    async def resubmit_with_more_resources(self, jg: bc.JobGroup) -> Optional[hb.Job]:
        global resources

        j = await jg._batch.get_job(self.job_id)
        status = await j.status()
        spec = status['spec']

        cpu = int(spec['resources']['req_cpu'])
        memory = spec['resources']['req_memory']

        resource_index = possible_resources.index((cpu, memory))
        if resource_index == -1 or resource_index == len(possible_resources) - 1:
            return None

        new_cpu, new_memory = possible_resources[resource_index + 1]

        resources = {
            'cpu': str(new_cpu),
            'memory': new_memory,
            'storage': spec['resources']['req_storage'],
            'preemptible': spec['resources']['preemptible'],
        }

        env = {var['name']: var['value'] for var in spec['env']}

        attributes = copy.deepcopy(status['attributes'])

        inputs = spec.get('input_files')
        if inputs:
            inputs = [(transfer['from'], transfer['to']) for transfer in inputs]

        outputs = spec.get('output_files')
        if outputs:
            outputs = [(transfer['from'], transfer['to']) for transfer in outputs]

        cloudfuse = spec.get('gcsfuse')
        if cloudfuse is not None:
            cloudfuse = [(mount['bucket'], mount['mount_path'], mount['read_only']) for mount in spec['gcsfuse']]

        j = jg.create_job(spec['process']['image'],
                          spec['process']['command'],
                          env=env,
                          resources=resources,
                          attributes=attributes,
                          input_files=inputs,
                          output_files=outputs,
                          always_run=status['always_run'],
                          timeout=spec.get('timeout'),
                          cloudfuse=cloudfuse,
                          requester_pays_project=spec.get('requester_pays_project'),
                          regions=spec['regions'],
                          always_copy_output=spec['always_copy_output'],
                          )

        return j


def heal(contig: str, billing_project: str, remote_tmpdir: str, max_attempts: int = 2):
    import os
    import re
    import asyncio
    from collections import Counter

    import hailtop.batch as hb

    batch_id = int(os.environ['HAIL_BATCH_ID'])
    job_group_id = int(os.environ['HAIL_JOB_GROUP_ID'])

    backend = hb.ServiceBackend(billing_project=billing_project, remote_tmpdir=remote_tmpdir)
    b = hb.Batch.from_batch_id(batch_id, backend=backend)
    jg = hb.JobGroup.from_job_group_id(b, job_group_id)

    PHASE_JOB_NAME_REGEX = re.compile(f'^phase/.*/{contig}/.*')

    async def _heal(b, jg):
        b_bc = b._async_batch
        jg_bc = jg._async_job_group

        while True:
            print('healing jobs...')

            job_info = [JobInfo.from_json(j) async for j in jg_bc.jobs()
                        if PHASE_JOB_NAME_REGEX.fullmatch(j['name']) is not None]

            print(f'found {len(job_info)} attempts')

            n_attempts = Counter(j.name for j in job_info)

            job_info.sort(key=lambda x: x.job_id, reverse=True)

            latest_attempts = {}
            for job in job_info:
                if job.name not in latest_attempts:
                    latest_attempts[job.name] = job
                    job.attempt_number = n_attempts[job.name]

            all_completed = all(j.state in ('Cancelled', 'Success') or (j.attempt_number >= max_attempts and j.state in ('Error', 'Failed'))
                                for j in latest_attempts.values())

            if all_completed:
                all_succeeded = all(j.state == 'Success' for j in latest_attempts.values())
                if all_succeeded:
                    return
                raise Exception('Some phasing jobs were cancelled, failed, or exceeded maximum number of attempts.')

            should_resubmit = False
            for j in latest_attempts.values():
                if j.attempt_number < max_attempts and j.should_be_rerun():
                    print(f'resubmitting job {j.job_id} {j.name} with attempt number {j.attempt_number + 1}')
                    await j.resubmit_with_more_resources(jg_bc)
                    should_resubmit |= True

            if should_resubmit:
                await b_bc.submit()

            await asyncio.sleep(60)

    asyncio.run(_heal(b, jg))


def heal_phase_jobs(b: hb.Batch,
                    jg: hb.JobGroup,
                    sample_group: SampleGroup,
                    contig: str,
                    docker: str,
                    billing_project: str,
                    remote_tmpdir: str,
                    max_attempts: int = 2) -> hb.Job:
    j = jg.new_python_job(name=f'phase-heal/sample-group-{sample_group.sample_group_index}/{contig}')
    j.image(docker)
    j.cpu(0.25)
    j.call(heal, contig, billing_project, remote_tmpdir, max_attempts)
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


def _vcf_to_mt(input_vcf: str, output_path: str):
    import hail as hl

    hl.init(backend="spark",
            local=f"local[{cpu}]",
            default_reference="GRCh38",
            tmp_dir="/io/",
            local_tmpdir="/io/",
            spark_conf={"spark.executor.memory": "7g", "spark.driver.memory": "7g", "spark.driver.maxResultSize": "7g"}
            )

    mt = hl.import_vcf(input_vcf)
    mt.write(output_path, overwrite=True)


def vcf_to_mt(b: hb.Batch,
              jg: hb.JobGroup,
              sample_group: SampleGroup,
              input_vcf: str,
              output_path: str,
              contig: str,
              docker: str,
              cpu: int,
              memory: str,
              storage: str,
              use_checkpoint: bool) -> Optional[Job]:
    if use_checkpoint and hfs.exists(output_path + '/_SUCCESS'):
        return None

    j = jg.new_python_job(attributes={'name': f'vcf-to-mt/sample-group-{sample_group.sample_group_index}/{contig}',
                                      'task': 'vcf-to-mt'})
    j.cpu(cpu)
    j.image(docker)
    j.storage(storage)
    j.memory(memory)

    j.call(_vcf_to_mt, input_vcf, output_path)

    return j


def _union(mt_paths: hb.ResourceFile,
           output_path: str,
           billing_project: str,
           remote_tmpdir: str,
           regions: str,
           contig: str,
           n_partitions: int):
    import subprocess as sp
    import hail as hl
    from hail.backend.service_backend import ServiceBackend

    setup = f'''
hailctl config set batch/billing_project "{billing_project}"
hailctl config set batch/remote_tmpdir "{remote_tmpdir}"
hailctl config set batch/regions "{','.join(regions)}"
'''

    sp.run(setup, capture_output=True, shell=True, check=True)

    hl.init(backend="batch", app_name=f"union-{contig}")

    backend = hl.current_backend()
    assert isinstance(backend, ServiceBackend)

    assert backend.regions == regions
    assert backend.remote_tmpdir == remote_tmpdir
    assert backend.billing_project == billing_project

    paths = []
    sample_sizes = []
    with open(mt_paths, 'r') as f:
        for line in f:
            path, sample_size = line.rstrip("\n").split('\t')
            paths.append(path)
            sample_sizes.append(int(sample_size))


    def add_info_if_needed(mt):
        return mt.annotate_rows(info=mt.info.annotate(INFO=mt.info.get("INFO", hl.null(hl.tarray(hl.tfloat64)))))


    mt_init = hl.read_matrix_table(paths[0])
    intervals = mt_init._calculate_new_partitions(n_partitions)

    mt_left = hl.read_matrix_table(paths[0], _intervals=intervals)
    mt_left = add_info_if_needed(mt_left)
    mt_left = mt_left.annotate_rows(
        info=mt_left.info.annotate(N=sample_sizes[0], AF=mt_left.info.AF[0], INFO=mt_left.info.INFO[0],
                                   RAF=mt_left.info.RAF[0]))
    mt_left = mt_left.annotate_rows(**{"info_0": mt_left.info})

    for idx, path in enumerate(paths[1:]):
        mt_right = hl.read_matrix_table(path, _intervals=intervals)
        mt_right = add_info_if_needed(mt_right)
        mt_right = mt_right.annotate_rows(
            info=mt_right.info.annotate(N=sample_sizes[idx], AF=mt_right.info.AF[0], INFO=mt_right.info.INFO[0],
                                        RAF=mt_right.info.RAF[0]))
        mt_left = mt_left.union_cols(mt_right,
                                     drop_right_row_fields=False,
                                     row_join_type='outer')

    mt = mt_left

    n_samples = mt.count_cols()
    n_batches = len(paths)

    mt = mt.annotate_rows(info=mt.info.annotate(AF=hl.array([mt[f"info_{i}"].AF for i in range(n_batches)])))
    mt = mt.annotate_rows(info=mt.info.annotate(INFO=hl.array([mt[f"info_{i}"].INFO for i in range(n_batches)])))
    mt = mt.annotate_rows(info=mt.info.annotate(N=hl.array([mt[f"info_{i}"].N for i in range(n_batches)])))


    def GLIMPSE_AF(mt):
        return hl.sum(hl.map(lambda af, n: af * n, mt.info.AF, mt.info.N)) / n_samples


    def GLIMPSE_INFO(mt):
        return hl.if_else((GLIMPSE_AF(mt) == 0) | (GLIMPSE_AF(mt) == 1),
                          1,
                          1 - hl.sum(
                              hl.map(lambda af, n, info: (1 - info) * 2 * n * af * (1 - af), mt.info.AF, mt.info.N,
                                     mt.info.INFO)) / (2 * n_samples * GLIMPSE_AF(mt) * (1 - GLIMPSE_AF(mt))))


    mt = mt.annotate_rows(info=mt.info.annotate(AF=GLIMPSE_AF(mt), INFO=GLIMPSE_INFO(mt)))

    mt = mt.annotate_rows(info=mt.info.drop('N'))
    mt = mt.drop(*[f'info_{i}' for i in range(n_batches)])
    mt = mt.drop(*[f'rsid_{i}' for i in range(1, n_batches)])
    mt = mt.drop(*[f'qual_{i}' for i in range(1, n_batches)])
    mt = mt.drop(*[f'filters_{i}' for i in range(1, n_batches)])

    if output_path.endswith('.mt'):
        mt.write(output_path)
        mt_count = hl.read_matrix_table(output_path)
        print(mt_count.count())
        print(mt.describe())
    else:
        assert output_path.endswith('.vcf.bgz')
        hl.export_vcf(mt, output_path, tabix=True)


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
                                  regions: List[str],
                                  use_checkpoints: bool,
                                  contig: str,
                                  n_partitions: int) -> Optional[Job]:
    if use_checkpoints:
        if output_path.endswith('.vcf.bgz') and hfs.exists(output_path):
            return None
        if output_path.endswith('.mt') and hfs.exists(output_path + '/_SUCCESS'):
            return None

    j = jg.new_python_job(attributes={'name': f'union/{contig}'})
    j.cpu(cpu)
    j.image(docker)
    j.storage(storage)
    j.memory(memory)
    j.spot(False)
    j.call(_union, vcf_paths, output_path, billing_project, remote_tmpdir, regions, contig, n_partitions)

    return j
