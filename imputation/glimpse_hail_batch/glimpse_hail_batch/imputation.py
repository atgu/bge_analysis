import argparse
import asyncio
from collections import defaultdict
import json
import re
from functools import partial

import hailtop.batch as hb
import hailtop.fs as hfs
from hailtop.aiotools.router_fs import RouterAsyncFS
from hailtop.utils import bounded_gather

from typing import Dict, List, Optional, Tuple

from .imputation_jobs import (
    copy_temp_crams_job, delete_temp_files_job,
    ligate, merge_vcfs, phase, union_sample_groups, write_success
)
from .globals import Chunk, SampleGroup, file_exists, find_crams, find_chunks, get_ligate_storage_requirement, split_samples_into_groups


async def run_sample_group(b: hb.Batch,
                           args: dict,
                           contig_chunks: Dict[str, List[Chunk]],
                           sample_group: SampleGroup,
                           fasta_input: hb.ResourceGroup,
                           ref_dict: hb.ResourceFile,
                           samples_per_copy_group: int,
                           prev_copy_cram_jobs: List[hb.Job],
                           fs: RouterAsyncFS) -> Tuple[List[hb.Job], Optional[hb.Job]]:
    print(f'staging sample group {sample_group.sample_group_index}')

    sample_group.write_sample_group_dict()

    skip_copy_crams = False
    skip_phasing = False
    skip_ligate = False
    skip_merge_vcf = False

    ligated_output_files_by_contig = sample_group.get_ligate_output_file_names(contig_chunks)
    phased_output_files = sample_group.get_phased_output_file_names(contig_chunks)

    phasing_already_completed = [False for contig, chunks in contig_chunks.items() for _ in chunks]

    if args['use_checkpoints']:
        phasing_already_completed = await bounded_gather(*[partial(file_exists, fs, phased_output_files[contig][chunk.chunk_idx] + '.bcf')
                                                                for contig, chunks in contig_chunks.items()
                                                                for chunk in chunks],
                                                                cancel_on_error=True)

        is_phasing_complete = all(phasing_already_completed)
        is_ligate_complete = all([hfs.exists(file + '.vcf.bgz') for contig, file in ligated_output_files_by_contig.items()])
        is_merge_complete = False # Check for success file hfs.is_dir(sample_group.merged_mt_path)

        skip_copy_crams = is_phasing_complete
        skip_phasing = is_phasing_complete
        skip_ligate = is_ligate_complete
        skip_merge_vcf = is_merge_complete

    copy_cram_jobs = []
    if not skip_copy_crams:
        for idx in range(0, len(sample_group.samples), samples_per_copy_group):
            # The CPU and memory here must be the same as the phase jobs in order
            # to avoid copying jobs happening much sooner than the subsequent phasing
            # jobs. By making the jobs the same size, they will be scheduled in subsequent
            # order.
            copy_j = copy_temp_crams_job(b,
                                         sample_group,
                                         idx,
                                         idx + samples_per_copy_group,
                                         args['phase_cpu'],
                                         args['phase_memory'])

            # this is really important to make sure that sample groups are copied when there's capacity for phasing jobs
            copy_j.depends_on(*prev_copy_cram_jobs)

            copy_cram_jobs.append(copy_j)

    phase_jobs = []

    n_variants_contig = {contig: sum(chunk.n_variants for chunk in chunks) for contig, chunks in contig_chunks.items()}

    if not skip_phasing:
        crams_list = sample_group.write_cram_list()
        crams_list_input = b.read_input(crams_list)

        phase_checkpoint_files = await bounded_gather(*[partial(sample_group.initialize_phased_glimpse_checkpoint_file, fs, contig, chunk.chunk_idx)
                                                        for contig, chunks in contig_chunks.items()
                                                        for chunk in chunks],
                                                      cancel_on_error=True)

        chunk_idx = 0
        for contig, chunks in contig_chunks.items():
            for chunk in chunks:
                phased_output_file = phased_output_files[contig][chunk.chunk_idx]

                phase_exists = phasing_already_completed[chunk_idx]
                phase_checkpoint_file = phase_checkpoint_files[chunk_idx]

                phase_j = phase(b,
                                phased_output_file,
                                phase_exists,
                                sample_group,
                                chunk,
                                sample_group.remote_cram_temp_dir,
                                sample_group.mount_path,
                                crams_list_input,
                                phase_checkpoint_file,
                                fasta_input,
                                args['docker_glimpse'],
                                args['phase_cpu'],
                                args['phase_memory'],
                                args['use_checkpoints'],
                                args['phase_impute_reference_only_variants'],
                                args['phase_call_indels'],
                                args['phase_n_burn_in'],
                                args['phase_n_main'],
                                args['phase_effective_population_size']
                                )

                if phase_j is not None:
                    phase_j.depends_on(*copy_cram_jobs)
                    phase_jobs.append(phase_j)

                chunk_idx += 1

    ligate_jobs = []
    if not skip_ligate:
        for contig, phased_files in phased_output_files.items():
            phased_inputs = [b.read_input_group(bcf=file + '.bcf', csi=file + '.bcf.csi')
                             for file in phased_output_files[contig]]

            ligate_storage_required = args['ligate_storage'] or get_ligate_storage_requirement(20, len(sample_group.samples), n_variants_contig[contig])

            ligate_j = ligate(b,
                              sample_group,
                              contig,
                              args['docker_glimpse'],
                              args['ligate_cpu'],
                              args['ligate_memory'],
                              ligate_storage_required,
                              phased_inputs,
                              ligated_output_files_by_contig[contig],
                              ref_dict,
                              args['use_checkpoints'])

            if ligate_j is not None:
                ligate_j.depends_on(*(copy_cram_jobs + phase_jobs))
                ligate_jobs.append(ligate_j)

    merge_vcf_jobs = []
    if not skip_merge_vcf:
        merge_j = merge_vcfs(b,
                             sample_group,
                             sample_group.ligate_output_dir + '*.vcf.bgz',
                             sample_group.merged_mt_path,
                             args['docker_hail'],
                             args['merge_vcf_cpu'],
                             args['merge_vcf_memory'],
                             args['merge_vcf_storage'],
                             args['use_checkpoints'])
        if merge_j is not None:
            merge_j.depends_on(*(copy_cram_jobs + phase_jobs + ligate_jobs))
            merge_vcf_jobs.append(merge_j)

    success_j = None
    if not args['use_checkpoints'] or not hfs.exists(sample_group.success_file):
        success_j = write_success(b, sample_group, args['docker_hail'])
        success_j.depends_on(*(copy_cram_jobs + phase_jobs + ligate_jobs + merge_vcf_jobs))

    if args['always_delete_temp_files'] or success_j is not None:
        delete_jobs = []
        delete_j = delete_temp_files_job(b, sample_group, args['save_checkpoints'])
        delete_j.depends_on(*(copy_cram_jobs + phase_jobs + ligate_jobs + merge_vcf_jobs))
        delete_j.always_run(True)
        delete_jobs.append(delete_j)

    return (copy_cram_jobs, success_j)


async def impute(args: dict):
    if hfs.exists(args['output_file']):
        raise Exception(f'output file {args["output_file"]} already exists.')

    with hfs.open(args['staging_remote_tmpdir'].rstrip('/') + '/config.json', 'w') as f:
        f.write(json.dumps(args, indent=4) + '\n')

    batch_regions = args['batch_regions']
    if batch_regions is not None:
        batch_regions = batch_regions.split(',')

    batch_name = args['batch_name'] or 'glimpse-imputation'

    backend = hb.ServiceBackend(billing_project=args['billing_project'],
                                remote_tmpdir=args['batch_remote_tmpdir'],
                                regions=batch_regions,
                                gcs_requester_pays_configuration=args['gcs_requester_pays_configuration'])

    if args['batch_id'] is not None:
        b = hb.Batch.from_batch_id(args['batch_id'], backend=backend, requester_pays_project=args['gcs_requester_pays_configuration'])
    else:
        b = hb.Batch(name=batch_name, backend=backend, requester_pays_project=args['gcs_requester_pays_configuration'])

    mount_point = '/crams/'

    samples = find_crams(args['sample_manifest'],
                         args['sample_id_col'],
                         args['cram_path_col'],
                         args['cram_index_path_col'],
                         args['n_samples'])

    sample_groups = split_samples_into_groups(samples,
                                              args['sample_group_size'],
                                              args['staging_remote_tmpdir'],
                                              mount_point)
    if args['sample_group_index'] is not None:
        sample_groups = [sg for sg in sample_groups if sg.sample_group_index == args['sample_group_index']]
        assert sample_groups

    chunks = find_chunks(args['reference_dir'],
                         args['chunk_info_dir'],
                         re.compile(args['binary_reference_file_regex']),
                         requested_contig=args['contig'],
                         requested_chunk_index=args['chunk_index'],
                         requester_pays_config=args['gcs_requester_pays_configuration'])

    chunks = [chunk for chunk in chunks if chunk.contig != "chrX"]

    print(f'found {len(chunks)} chunks')
    print(f'found {len(sample_groups)} sample groups')

    contig_chunks = defaultdict(list)
    for chunk in chunks:
        contig_chunks[chunk.contig].append(chunk)

    fasta_input = b.read_input_group(**{'fasta': args['fasta'], 'fasta.fai': f'{args["fasta"]}.fai'})
    ref_dict = b.read_input(args['ligate_ref_dict'])

    success_jobs = []
    prev_copy_cram_jobs = []
    for sample_group in sample_groups:
        prev_copy_cram_jobs, success_j = await run_sample_group(b,
                                                                args,
                                                                contig_chunks,
                                                                sample_group,
                                                                fasta_input,
                                                                ref_dict,
                                                                args['samples_per_copy_group'],
                                                                prev_copy_cram_jobs,
                                                                backend._fs)
        if success_j is not None:
            success_jobs.append(success_j)

    sample_group_merged_mts = [sample_group.merged_mt_path for sample_group in sample_groups]
    union_sample_groups_inputs_path = args['staging_remote_tmpdir'].rstrip('/') + '/sample_group_mts.txt'
    with hfs.open(union_sample_groups_inputs_path, 'w') as f:
        for merged_mt in sample_group_merged_mts:
            f.write(merged_mt + '\n')

    union_j = union_sample_groups(b,
                                  b.read_input(union_sample_groups_inputs_path),
                                  args['output_file'],
                                  args['docker_hail'],
                                  2,
                                  'standard',
                                  '10Gi',
                                  batch_name,
                                  args['billing_project'],
                                  args['batch_remote_tmpdir'],
                                  batch_regions,
                                  args['use_checkpoints'],
                                  n_partitions=len(chunks))
    if union_j is not None:
        union_j.depends_on(*success_jobs)

    b.run(disable_progress_bar=False)

    backend.close()


if __name__ == '__main__':
    parser = argparse.ArgumentParser()

    parser.add_argument('--billing-project', type=str, required=False)
    parser.add_argument('--batch-remote-tmpdir', type=str, required=False)
    parser.add_argument('--batch-regions', type=str, default="us-central1")
    parser.add_argument('--batch-name', type=str, default='glimpse2-imputation')
    parser.add_argument('--batch-id', type=int, required=False)

    parser.add_argument('--seed', type=int, required=False, default=2423413432)

    parser.add_argument('--docker-glimpse', type=str, required=True)
    parser.add_argument('--docker-hail', type=str, required=True)

    parser.add_argument("--reference-dir", type=str, required=True)
    parser.add_argument("--chunk-info-dir", type=str, required=True)
    parser.add_argument('--binary-reference-file-regex', type=str, required=True)

    parser.add_argument('--sample-manifest', type=str, required=True)
    parser.add_argument('--sample-id-col', type=str, required=True)
    parser.add_argument('--cram-path-col', type=str, required=True)
    parser.add_argument('--cram-index-path-col', type=str, required=True)

    parser.add_argument('--n-samples', type=int, required=False)

    parser.add_argument('--sample-group-size', type=int, required=True, default=100)
    parser.add_argument('--samples-per-copy-group', type=int, default=100)

    parser.add_argument('--ligate-cpu', type=int, required=True)
    parser.add_argument('--ligate-memory', type=str, required=False, default='standard')
    parser.add_argument('--ligate-storage', type=str, required=False)
    parser.add_argument('--ligate-ref-dict', type=str, required=False, default="gs://gcp-public-data--broad-references/hg38/v0/Homo_sapiens_assembly38.dict")

    # Filters for running a specific set of samples and chunk
    parser.add_argument('--contig', type=str, required=False)
    parser.add_argument('--chunk-index', type=int, required=False)
    parser.add_argument('--sample-group-index', type=int, required=False)

    parser.add_argument('--staging-remote-tmpdir', type=str, required=True)
    parser.add_argument('--output-file', type=str, required=True)

    parser.add_argument('--fasta', type=str, required=False, default='gs://gcp-public-data--broad-references/hg38/v0/Homo_sapiens_assembly38.fasta')

    parser.add_argument('--use-checkpoints', action='store_true', required=False)
    parser.add_argument('--save-checkpoints', action='store_true', required=False)
    parser.add_argument('--always-delete-temp-files', action='store_true', required=False)

    # Extra phase arguments
    parser.add_argument('--phase-cpu', type=int, required=True)
    parser.add_argument('--phase-memory', type=str, required=False, default='standard')
    parser.add_argument('--phase-impute-reference-only-variants', action='store_true', required=False)
    parser.add_argument('--phase-call-indels', action='store_true', required=False)
    parser.add_argument('--phase-n-burn-in', type=int, required=False)
    parser.add_argument('--phase-n-main', type=int, required=False)
    parser.add_argument('--phase-effective-population-size', type=int, required=False)

    # Extra merge vcf arguments
    parser.add_argument('--merge-vcf-cpu', type=int, required=True)
    parser.add_argument('--merge-vcf-memory', type=str, required=False, default='lowmem')
    parser.add_argument('--merge-vcf-storage', type=str, required=False, default='0Gi')

    parser.add_argument('--gcs-requester-pays-configuration', type=str, required=False)

    args = vars(parser.parse_args())

    print('submitting jobs with the following parameters:')
    print(json.dumps(args, indent=4))
    asyncio.run(impute(args))
