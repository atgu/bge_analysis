import argparse
import os
import re
from collections import defaultdict
from typing import List

import hailtop.batch as hb
import hailtop.fs as hfs


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
                                regions=batch_regions,
                                gcs_requester_pays_configuration=args['gcs_requester_pays_configuration'])

    b = hb.Batch(name=batch_name, backend=backend, requester_pays_project=args['gcs_requester_pays_configuration'])

    reference_dir = args['reference_dir']
    chunk_info_dir = args['chunk_info_dir']
    binary_reference_file_regex = re.compile(args['binary_reference_file_regex'])
    chunk_file_regex = re.compile(args['chunk_file_regex'])

    reference_files = hfs.ls(reference_dir, requester_pays_config=args['gcs_requester_pays_configuration'])

    if chunk_info_dir is None:
        contig_chunk_info = defaultdict(list)
        for reference_file in reference_files:
            match = binary_reference_file_regex.fullmatch(os.path.basename(reference_file.path)).groupdict()

            chunk_idx = int(match['chunk_index'])
            contig = match['contig']
            start = int(match['start'])
            end = int(match['end'])

            contig_chunk_info[contig].append((chunk_idx, start, end))

        chunk_info = {}
        for contig, chunks in contig_chunk_info.items():
            dummy_chunk_path = args['batch_remote_tmpdir'] + f'/temp_chunk_files/{contig}_chunks.txt'
            with hfs.open(dummy_chunk_path, 'w') as f:
                for chunk_idx, start, end in chunks:
                    data = [chunk_idx, contig, f'{contig}:{start}-{end}', f'{contig}:{start}-{end}', 0, 0, 0, 0]
                    data = "\t".join([str(d) for d in data])
                    f.write(f'{data}\n')
            chunk_info[contig] = dummy_chunk_path
    else:
        chunk_info_files = hfs.ls(chunk_info_dir, requester_pays_config=args['gcs_requester_pays_configuration'])

        chunk_info = {chunk_file_regex.fullmatch(os.path.basename(chunk_info_file.path)).groupdict()['contig']: chunk_info_file.path
                      for chunk_info_file in chunk_info_files}

    contig_chunk_info = defaultdict(list)

    for file in reference_files:
        reference_chunk = b.read_input(file.path)
        match = binary_reference_file_regex.fullmatch(os.path.basename(file.path)).groupdict()

        chunk_index = int(match['chunk_index'])
        contig = match['contig']

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
    parser.add_argument('--gcs-requester-pays-configuration', type=str, required=False)

    parser.add_argument('--docker-glimpse-extract-sites', type=str, required=False,
                        default="us.gcr.io/broad-dsde-methods/glimpse_extract_num_sites_from_reference_chunks:michaelgatzen_edc7f3a")

    parser.add_argument("--reference-dir", type=str, required=True)
    parser.add_argument("--chunk-info-dir", type=str, required=False)
    parser.add_argument('--binary-reference-file-regex', type=str, required=True)

    parser.add_argument('--output-directory', type=str, required=True)

    args = vars(parser.parse_args())
    regenerate_chunk_metadata(args)
