import argparse
from collections import defaultdict
from typing import List

import hailtop.batch as hb
import hailtop.fs as hfs

from .globals import chunk_info_file_dir_str, parse_binary_reference_path, parse_chunk_path, reference_file_dir_str


def get_number_of_sites_in_chunk(b: hb.Batch,
                                 reference_chunk: hb.ResourceFile,
                                 docker: str,
                                 contig: str,
                                 chunk_idx: int,
                                 original_chunk_info: hb.ResourceFile):
    j = b.new_bash_job(name='get-number-of-sites')
    j.image(docker)
    j.cpu(1)
    j.memory('standard')

    cmd = f'''
/bin/GLIMPSE2_extract_num_sites_from_reference_chunk {reference_chunk} > n_sites.txt
cat n_sites.txt
N_RARE=`grep "Lrare" n_sites.txt | sed 's/Lrare=//'`
N_COMMON=`grep "Lcommon" n_sites.txt | sed 's/Lcommon=//'`
cat {original_chunk_info} | awk '{{if ($1 == "{chunk_idx}" && $2 == "{contig.split('_')[0]}") print $0 }}' | \
    awk -v n_rare=$N_RARE -v n_common=$N_COMMON -F'\t' 'BEGIN {{OFS = FS}} {{$7 = n_rare + n_common; $8 = n_common; print}}' > {j.new_chunk_info}
cat {j.new_chunk_info}
'''

    j.command(cmd)

    return j


def concatenate(b: hb.Batch, contig: str, files: List[hb.ResourceFile]):
    j = b.new_bash_job(name=f'concat-{contig}')
    j.cpu(1)
    j.memory('standard')
    j.command('cat ' + ' '.join(files) + f' | sort -nk 1 > {j.chunk_info}')
    return j


def regenerate_chunk_metadata(args: dict):
    batch_regions = args['batch_regions']
    if batch_regions is not None:
        batch_regions = batch_regions.split(',')

    batch_name = args['batch_name'] or 'regenerate-chunk-metadata'

    backend = hb.ServiceBackend(billing_project=args['billing_project'],
                                remote_tmpdir=args['batch_remote_tmpdir'],
                                regions=batch_regions)
    b = hb.Batch(name=batch_name, backend=backend)

    reference_dir = reference_file_dir_str(args['split_reference_dir'])
    reference_files = hfs.ls(reference_dir)

    chunk_info_dir = chunk_info_file_dir_str(args['split_reference_dir'])
    chunk_info_files = hfs.ls(chunk_info_dir)

    chunk_info = {parse_chunk_path(chunk_info_file.path): chunk_info_file.path
                  for chunk_info_file in chunk_info_files}

    contig_chunk_info = defaultdict(list)

    for file in reference_files:
        reference_chunk = b.read_input(file.path)
        contig, chunk_index = parse_binary_reference_path(file.path)

        original_chunk_info = b.read_input(chunk_info[contig])

        j = get_number_of_sites_in_chunk(b,
                                         reference_chunk,
                                         args['docker_glimpse_extract_sites'],
                                         contig,
                                         chunk_index,
                                         original_chunk_info)

        contig_chunk_info[contig].append(j.new_chunk_info)

    for contig, new_chunk_info_files in contig_chunk_info.items():
        j = concatenate(b, contig, new_chunk_info_files)
        b.write_output(j.chunk_info, args['output_directory'].rstrip('/') + f'/chunks_{contig}.txt')

    b.run()


if __name__ == '__main__':
    parser = argparse.ArgumentParser()

    parser.add_argument('--billing-project', type=str, required=False)
    parser.add_argument('--batch-remote-tmpdir', type=str, required=False)
    parser.add_argument('--batch-regions', type=str, default="us-central1")
    parser.add_argument('--batch-name', type=str, default='regenerate-chunk-info')

    parser.add_argument('--docker-glimpse-extract-sites', type=str, required=False,
                        default="us.gcr.io/broad-dsde-methods/glimpse_extract_num_sites_from_reference_chunks:michaelgatzen_edc7f3a")

    parser.add_argument("--split-reference-dir", type=str, required=True)

    parser.add_argument('--output-directory', type=str, required=True)

    args = vars(parser.parse_args())
    regenerate_chunk_metadata(args)
