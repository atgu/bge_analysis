import argparse
import hailtop.batch as hb
import hailtop.fs as hfs
from typing import List, Tuple, Optional
import math
import json
import subprocess as sp


from .globals import chunks_file_str, chunk_info_file_dir_str, reference_file_dir_str, reference_file_str


def file_size(file: str) -> int:
    result = sp.run(f"gcloud storage ls {file} --json", shell=True, check=True, capture_output=True).stdout
    result = json.loads(result)
    return int(result[0]["metadata"]["size"])


def read_regions_file(regions_file: str, requested_region: Optional[str]) -> Tuple[List[str], List[str], List[str], List[str]]:
    contig_regions = []
    reference_panel_filenames = []
    genetic_map_filenames = []

    with hfs.open(regions_file, 'r') as f:
        f.readline()
        for line in f:
            line = str(line)
            line = line.rstrip()
            region, reference_filename, genetic_map_filename = line.split('\t')
            contig_regions.append(region)
            reference_panel_filenames.append(reference_filename)
            genetic_map_filenames.append(genetic_map_filename)

    if requested_region is not None:
        trunc_region = requested_region.split(':')[0]
        idx = contig_regions.index(trunc_region)
        if idx != -1:
            contig_regions = [requested_region]
            reference_panel_filenames = [reference_panel_filenames[idx]]
            genetic_map_filenames = [genetic_map_filenames[idx]]
        else:
            raise Exception(f"unknown region: {region}")

    reference_panel_sizes = []
    for file in reference_panel_filenames:
        size = file_size(file)
        reference_panel_sizes.append(size)

    genetic_map_sizes = []
    for file in genetic_map_filenames:
        size = file_size(file)
        genetic_map_sizes.append(size)

    region_storage_sizes = [2.2 * ref_size + gm_size + 10 for ref_size, gm_size in zip(reference_panel_sizes, genetic_map_sizes)]
    region_storage_sizes = [f'{math.ceil(size / 1024 / 1024 / 1024)}Gi' for size in region_storage_sizes]

    return (contig_regions, reference_panel_filenames, genetic_map_filenames, region_storage_sizes)


def split_reference_job(b: hb.Batch,
                        args: dict,
                        contig_index: int,
                        contig: str,
                        reference_panel: hb.ResourceGroup,
                        genetic_map: hb.ResourceFile,
                        region_storage_size: str):
    split_reference_output_dir = args['output_dir'].rstrip('/')

    split_reference_dir = reference_file_dir_str(split_reference_output_dir)
    chunks_dir = chunk_info_file_dir_str(split_reference_output_dir)
    chunks_file = chunks_file_str(chunks_dir, contig)

    j = b.new_job(name=f'split-reference-{contig_index}-{contig}',
                  attributes={'contig-index': str(contig_index),
                              'contig': str(contig)})
    j.image(args['docker'])
    j.cpu(args['cpu'])
    j.memory(args['memory'])
    j.storage(region_storage_size)

    glimpse_chunk_cmd = f"""
/bin/GLIMPSE2_chunk \
    --input {reference_panel["bcf"]} \
    --region {contig} \
    --map {genetic_map} \
    --sequential \
    --threads {args['cpu']} \
    --output ${{CHUNKS_FILE}}"""

    if args['seed'] is not None:
        glimpse_chunk_cmd += f' \ --seed {args["seed"]}'
    if args['min_window_cm'] is not None:
        glimpse_chunk_cmd += f' \ --window-cm {args["min_window_cm"]}'
    if args['uniform_number_of_variants'] is not None:
        glimpse_chunk_cmd += f' \ --uniform-number-variants'

    glimpse_split_ref_cmd = f'''
/bin/GLIMPSE2_split_reference \
    --threads {args['cpu']} \
    --reference {reference_panel["bcf"]} \
    --map {genetic_map} \
    --input-region ${{IRG}} \
    --output-region ${{ORG}} \
    --output ${{REFERENCE_OUTPUT_DIR}}/{reference_file_str()}'''

    if args['keep_monomorphic_ref_sites']:
        glimpse_split_ref_cmd += f' \ --keep-monomorphic-ref-sites'
    if args['seed'] is not None:
        glimpse_split_ref_cmd += f' \ --seed {args["seed"]}'

    j.command(f'''
set -xeuo pipefail

# Print chunk index to variable
CONTIGINDEX=$(printf "%04d" {contig_index})
CONTIG="{contig}"

CHUNKS_FILE="{j.chunks}"
REFERENCE_OUTPUT_DIR="{j.reference_output_dir}"

# this is a hack to make sure the index is newer than the bcf file
touch {reference_panel["bcf.csi"]}

{glimpse_chunk_cmd}

if [ -f chunks_contigindex_${{CONTIGINDEX}}.txt_uniform ]; then
    mv chunks_contigindex_${{CONTIGINDEX}}.txt_uniform ${{CHUNKS_FILE}}
fi

mkdir -p ${{REFERENCE_OUTPUT_DIR}}

I_CHUNK=0
while IFS="" read -r LINE || [ -n "$LINE" ];
do
    # Extract coordinates from chunks.txt file
    printf -v ID "%02d" $(echo $LINE | cut -d" " -f1)
    IRG=$(echo $LINE | cut -d" " -f3)
    ORG=$(echo $LINE | cut -d" " -f4)

    # Print chunk index to variable
    CHUNKINDEX=$(printf "%04d" $I_CHUNK)

    {glimpse_split_ref_cmd}

    # Increase i (and make sure the exit code is zero)
    (( I_CHUNK++ )) || true
done < ${{CHUNKS_FILE}}
''')

    b.write_output(j.chunks, chunks_file)
    b.write_output(j.reference_output_dir, split_reference_dir)

    return j


def split_reference(backend: hb.ServiceBackend, args: dict):
    b = hb.Batch('glimpse-split-reference', backend=backend)

    regions_file = args.pop('regions_file')
    contig_regions, reference_panel_filenames, genetic_map_filenames, region_storage_sizes = read_regions_file(regions_file, args["region"])

    for contig_idx, (contig, reference_panel, genetic_map, region_storage) in enumerate(zip(contig_regions, reference_panel_filenames, genetic_map_filenames, region_storage_sizes)):
        reference_panel = b.read_input_group(**{"bcf": reference_panel, "bcf.csi": f'{reference_panel}.csi'})
        genetic_map = b.read_input(genetic_map)
        split_reference_job(b, args, contig_idx, contig, reference_panel, genetic_map, region_storage)

    b.run()

if __name__ == '__main__':
    p = argparse.ArgumentParser()

    p.add_argument('--billing-project', type=str, required=False)
    p.add_argument('--remote-tmpdir', type=str, required=False)
    p.add_argument('--batch-regions', type=str, default="us-central1")

    p.add_argument("--regions-file", type=str, required=True)
    p.add_argument("--region", type=str, required=False)

    p.add_argument('--seed', type=int, default=3242342)
    p.add_argument("--min-window-cm", type=int, required=False)
    p.add_argument("--uniform-number-of-variants", action="store_true", default=False)
    p.add_argument("--keep-monomorphic-ref-sites", action="store_true", default=False)

    p.add_argument("--docker", type=str, required=True)
    p.add_argument('--cpu', type=int, required=False, default=4)
    p.add_argument('--memory', type=str, required=False, default="standard")

    p.add_argument('--output-dir', type=str, required=True)
    p.add_argument('--batch-name', type=str, required=False)

    args = p.parse_args()

    batch_regions = args.batch_regions
    if batch_regions:
        batch_regions = args.batch_regions.split(',')

    backend = hb.ServiceBackend(billing_project=args.billing_project,
                                remote_tmpdir=args.remote_tmpdir,
                                regions=batch_regions)

    split_reference(backend, vars(args))

    backend.close()
