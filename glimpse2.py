__author__ = 'Mary T. Yohannes'

### GLIMPSE2 SCRIPT #1 ###

# In this script, you'll find how to (using GLIMPSE2 and bcftools):
# # Convert the reference panel into GLIMPSE2’s binary file format
# # Impute and phase a whole chromosome
# # Ligate imputed chunks of the the same chromosome
### Based on the tutorial probided here https://odelaneau.github.io/GLIMPSE/docs/tutorials/getting_started/

import hailtop.batch as hb
import subprocess
import os
import pandas as pd
import hail as hl


# get file size to allocate resources for batch jobs accordingly
def get_file_size(file):
    file_info = hl.utils.hadoop_stat(file)
    size_bytes = file_info['size_bytes']
    size_gigs = size_bytes / (1024 * 1024 * 1024)
    return size_gigs

# Step 1) Convert the reference panel into GLIMPSE2’s binary file format (one per chunk)
def create_binary_ref(ref_bcf, genetic_map, chr, regs, rege, irg, org):
    j = b.new_job(name=f'create binary ref panel - {chr}:{regs}-{rege}')  # define job
    j.cpu(4)
    j.storage('10Gi')
    j.image('docker.io/simrub/glimpse:v2.0.0-27-g0919952_20221207')  # Docker image with GLIMPSE2

    j.command(f''' GLIMPSE2_split_reference --reference {ref_bcf['bcf']} \
        --map {genetic_map} \
        --input-region {irg} \
        --output-region {org} \
        --threads 4 \
        --output binary_reference_panel 
        
mv binary_reference_panel_{chr}_{regs}_{rege}.bin {j.ofile}''')
    return j

# Step 2) Impute and phase the chunks that are in a binary format
def glimpse2_phase(cram_list, ref_binary, chr, regs, rege, storage_size):
    j = b.new_job(name=f'run GLIMPSE2_phase - {chr}:{regs}-{rege}')  # define job
    j.cpu(4)
    j.storage(f'{storage_size}Gi')  # increase storage
    j.memory('highmem')
    j.image('docker.io/simrub/glimpse:v2.0.0-27-g0919952_20221207')  # Docker image with GLIMPSE2 and bcftools

    # output handler
    j.declare_resource_group(ofile={
        'bcf': '{root}.bcf',
        'bcf.csi': '{root}.bcf.csi'})

    # the input file is a list of cram files so get the names/paths of those files
    cram_txt = '\n'.join([f'{c["cram"]}' for c in cram_list])

    # create a text file containing the full list of bcf files, one imputed chunk file per line
    j.command(f''' echo "{cram_txt}" > cram_list.txt''')

    # run GLIMPSE2_phase
    j.command(f''' GLIMPSE2_phase --bam-list cram_list.txt \
    --reference {ref_binary} \
    --threads 4 \
    --output {j.ofile['bcf']} ''')

    # generate index file
    j.command(f''' bcftools index -f {j.ofile['bcf']} ''')
    return j

# Step 3) Ligate imputed chunks of the the same chromosome
def ligate_chunks(imputed_chunk_list, chr):
    j = b.new_job(name=f'ligate chunks - {chr}')  # define job
    j.cpu(2)
    j.image('docker.io/simrub/glimpse:v2.0.0-27-g0919952_20221207')  # Docker image with GLIMPSE and bcftools

    # output handler
    j.declare_resource_group(ofile={
        'bcf': '{root}.bcf',
        'bcf.csi': '{root}.bcf.csi'})

    # the input file is a list of imputed chunks - have to be in order
    # here we are getting the names/paths of those chunks
    imp_chunk_bcf_names = '\n'.join([f'{c["bcf"]}' for c in imputed_chunk_list])

    # create a text file containing the full list of bcf files, one imputed chunk file per line
    j.command(f''' echo "{imp_chunk_bcf_names}" > chunk_list.txt''')

    j.command(f''' GLIMPSE2_ligate --input chunk_list.txt --output {j.ofile['bcf']} --threads 2 ''')

    # generate index file
    j.command(f''' bcftools index -f {j.ofile['bcf']} ''')
    return j


if __name__ == '__main__':
    backend = hb.ServiceBackend(billing_project='diverse-pop-seq-ref', remote_tmpdir='gs://bge/glimpse2') # set up backend

    # do this for each chromosome - one batch per chr
    for n in range(1, 23):

        b = hb.Batch(backend=backend, name=f'GLIMPSE2-chr{n}-impute&ligate-20KEMRI')  # define batch for each chromosome

        # a reference panel of haplotypes (only bi-allelic snps) -  files obtained from running ref_qc.py
        bcf_path = f'gs://bge/reference_panel/hgdp.tgp.trio/hgdp.tgp.chr{n}.onlysnps.bcf'
        ref_bcf = b.read_input_group(**{'bcf': bcf_path,
                                        'bcf.csi': f'{bcf_path}.csi'})

        # a fine-scale genetic map - b38
        genetic_map = b.read_input(f'gs://hgdp-1kg/phasing/maps/b38/chr{n}.b38.gmap.gz')

        # file containing the imputation chunks per chr - obtained from GLIMPSE developers (4cM)
        chunk_file = pd.read_csv(f'gs://hgdp-1kg/phasing/chunks/b38/4cM/chunks_chr{n}.txt', sep='\t', header=None)
        regions = [(irg, org) for irg, org in zip(chunk_file[2], chunk_file[3])]  # grab buffered and imputation regions from txt file

        results = []  # empty list for imputed chunks - STEP 2 glimpse2_phase outputs

        # for each buffered and imputation regions, do the following
        for i in range(len(regions)):

            regs = regions[i][0].split(':')[1].split('-')[0] # start region
            rege = regions[i][0].split(':')[1].split('-')[1] # end region

            # RUN STEP 1 - create binary reference panel
            run_binary_ref = create_binary_ref(ref_bcf, genetic_map, f'chr{n}', regs, rege, regions[i][0], regions[i][1])

            # get cram file google cloud paths for target samples
            root_path = 'gs://fc-d6e8ddd0-6b00-4f7c-8051-351ab37f79c1/NeuroGAP-Psychosis_KEMRI_Psychosis_BGE/RP-1889.PDO-29432/Exome'  # common root path
            sub_path = '**/**/*.cram'
            cram_paths = subprocess.check_output(['gsutil', 'ls', f'{root_path}/{sub_path}']).decode().strip('\n').split('\n')  # returns a list of individual cram file paths

            # USE ONLY 20 SAMPLES FOR NOW
            cram_paths = cram_paths[0:20]

            cfiles = [] # list to hold cram files - imput for STEP 2 glimpse2_phase function
            total_size = 0
            for path in cram_paths:
                # only file name - without path (removed by os.path.basename) and .cram (removed by os.path.splitext) - for keeping track
                #label = os.path.splitext(os.path.basename(path))[0]  # ex output: 15001708095242PK

                # read in both the .cram and .crai files
                sample_cram = b.read_input_group(
                    cram=path,
                    crai=f'{path}.crai')

                cfiles.append(sample_cram) # list of cram file paths

                # for job storage size
                bcf_size = round(get_file_size(path))
                total_size = total_size + bcf_size

            # RUN STEP 2  - Running GLIMPSE2
            run_glimpse2_phase = glimpse2_phase(cfiles, run_binary_ref.ofile, f'chr{n}', regs, rege, total_size)

            results.append(run_glimpse2_phase.ofile)  # append outputs to the list 'results'

        # RUN STEP 3 - ligate chunks of the the same chromosome
        run_ligate_chunks = ligate_chunks(results, f'chr{n}')

        # write out merged GLs for each chr
        b.write_output(run_ligate_chunks.ofile, f'gs://bge/glimpse2/20KEMRI/glimpse2_imputation/imputed.ligated.chr{n}')

        b.run(wait=False) # run batch

    backend.close()   # close a batch backend







