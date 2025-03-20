import argparse
import asyncio
import base64
import json
import os
from shlex import quote as shq

import typer
from rich.prompt import Confirm, IntPrompt, Prompt
from typer import Abort, Exit

from hailtop.aiotools.copy import copy_from_dict
import hailtop.batch as hb
from hailtop.utils import secret_alnum_string


async def submit(args):
    tmpdir_path_prefix = secret_alnum_string()

    remote_tmpdir = args['batch_remote_tmpdir']
    billing_project = args['billing_project']
    regions = args['batch_regions'].split(',')

    def cloud_prefix(path):
        return f'{remote_tmpdir}/{tmpdir_path_prefix}/{os.path.basename(path)}'

    sample_manifest_cloud_file = cloud_prefix(args['sample_manifest'])

    backend = hb.ServiceBackend(billing_project=billing_project, regions=regions, remote_tmpdir=remote_tmpdir)

    if args['batch_id'] is not None:
        b = hb.Batch.from_batch_id(args['batch_id'], backend=backend)
    else:
        b = hb.Batch(name=args['batch_name'], backend=backend)

    j = b.new_bash_job(name='submit-jobs')
    j.image(args['docker_hail'])

    await copy_from_dict(
        files=[
            {'from': args['sample_manifest'], 'to': sample_manifest_cloud_file}
        ]
    )

    local_sample_manifest = '/sample_manifest.tsv'
    args['sample_manifest'] = local_sample_manifest

    arguments_str = base64.b64encode(json.dumps(args).encode('utf-8')).decode('utf-8')

    sample_manifest_input = b.read_input(sample_manifest_cloud_file)

    j.command(f'mv {sample_manifest_input} {local_sample_manifest}')
    j.command(f'python3 -m glimpse_hail_batch.imputation.imputation "{shq(arguments_str)}"')

    batch_handle = await b._async_run(wait=False, disable_progress_bar=True)
    assert batch_handle
    print(f'Submitted batch {batch_handle.id}')

    await backend.async_close()


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
    parser.add_argument('--chunk-file-regex', type=str, required=True)

    parser.add_argument('--sample-manifest', type=str, required=True)
    parser.add_argument('--sample-id-col', type=str, required=True)
    parser.add_argument('--cram-path-col', type=str, required=True)
    parser.add_argument('--cram-index-path-col', type=str, required=True)
    parser.add_argument('--sex-col', type=str, required=False)
    parser.add_argument('--female-code', type=str, required=False)

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
    parser.add_argument('--phase-memory', type=str, required=True)
    parser.add_argument('--phase-impute-reference-only-variants', action='store_true', required=False)
    parser.add_argument('--phase-call-indels', action='store_true', required=False)
    parser.add_argument('--phase-n-burn-in', type=int, required=False)
    parser.add_argument('--phase-n-main', type=int, required=False)
    parser.add_argument('--phase-effective-population-size', type=int, required=False)
    parser.add_argument('--phase-max-attempts', type=int, required=False, default=2)

    # Extra merge vcf arguments
    parser.add_argument('--merge-vcf-cpu', type=int, required=True)
    parser.add_argument('--merge-vcf-memory', type=str, required=False, default='lowmem')
    parser.add_argument('--merge-vcf-storage', type=str, required=False, default='0Gi')

    parser.add_argument('--gcs-requester-pays-configuration', type=str, required=False)

    parser.add_argument('--non-par-contigs', type=str, required=False)

    args = vars(parser.parse_args())

    print('submitting jobs with the following parameters:')
    print(json.dumps(args, indent=4))

    if args['phase_memory'] == 'highmem':
        correct_memory = Confirm.ask('Are you sure you want to use highmem machines for phasing?')
        if not correct_memory:
            raise Exit()

    asyncio.run(submit(args))
