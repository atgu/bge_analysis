__author__ = 'Mary T. Yohannes'

### GLIMPSE2 SCRIPT #1 ###

# In this script, you'll find how to (using GLIMPSE2 and bcftools):
# # Convert the reference panel into GLIMPSE2’s binary file format
# # Impute and phase chromosome chunks using the genotype likelihoods generated for GLIMPSE1 and --ne 20000
# # Ligate the imputed chunks back together per chromosome
### Based on the tutorial provided here (a little change - used --input-gl option instead of --bam-list for step 2) https://odelaneau.github.io/GLIMPSE/docs/tutorials/getting_started/

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
    j.image(
        'docker.io/imary116/glimpse2:with-bcftools-and-updated-info-score')  # Docker image with GLIMPSE2 with INFO score corrected
    # j.image('docker.io/simrub/glimpse:v2.0.0-27-g0919952_20221207')  # Docker image with GLIMPSE2

    j.command(f''' GLIMPSE2_split_reference --reference {ref_bcf['bcf']} \
        --map {genetic_map} \
        --input-region {irg} \
        --output-region {org} \
        --threads 4 \
        --output binary_reference_panel 

mv binary_reference_panel_{chr}_{regs}_{rege}.bin {j.ofile}''')
    return j


# Step 2) Impute and phase the chunks that are in a binary format
# use the genotype likelihoods that were generated for GLIMPSE1 - using --input-gl option instead of --bam-list
def glimpse2_phase(gl_vcf, ref_binary, chr, regs, rege, storage_size):
    j = b.new_job(name=f'run GLIMPSE2_phase - {chr}:{regs}-{rege}')  # define job
    j.cpu(4)
    j.storage(f'{storage_size}Gi')  # increase storage
    j.memory('highmem')
    j.image(
        'docker.io/imary116/glimpse2:with-bcftools-and-updated-info-score')  # Docker image with GLIMPSE2 with INFO score corrected
    # j.image('docker.io/simrub/glimpse:v2.0.0-27-g0919952_20221207')  # Docker image with GLIMPSE2 and bcftools

    # output handler
    j.declare_resource_group(ofile={
        'bcf': '{root}.bcf',
        'bcf.csi': '{root}.bcf.csi'})

    # run GLIMPSE2_phase
    j.command(f''' GLIMPSE2_phase --input-gl {gl_vcf['vcf.gz']} \
    --reference {ref_binary} \
    --ne 20000 \
    --threads 4 \
    --output {j.ofile['bcf']} ''')

    # generate index file
    j.command(f''' bcftools index -f {j.ofile['bcf']} ''')
    return j


# Step 3) Ligate imputed chunks of the the same chromosome
def ligate_chunks(imputed_chunk_list, chr, storage_size):
    j = b.new_job(name=f'ligate chunks - {chr}')  # define job
    j.cpu(4)
    j.storage(f'{storage_size}Gi')
    j.image(
        'docker.io/imary116/glimpse2:with-bcftools-and-updated-info-score')  # Docker image with GLIMPSE2 with INFO score corrected
    # j.image('docker.io/simrub/glimpse:v2.0.0-27-g0919952_20221207')  # Docker image with GLIMPSE and bcftools

    # output handler
    j.declare_resource_group(ofile={
        'bcf': '{root}.bcf',
        'bcf.csi': '{root}.bcf.csi'})

    # the input file is a list of imputed chunks - have to be in order
    # here we are getting the names/paths of those chunks
    imp_chunk_bcf_names = '\n'.join([f'{c["bcf"]}' for c in imputed_chunk_list])

    # create a text file containing the full list of bcf files, one imputed chunk file per line
    j.command(f''' echo "{imp_chunk_bcf_names}" > chunk_list.txt''')

    j.command(f''' GLIMPSE2_ligate --input chunk_list.txt --output {j.ofile['bcf']} --threads 4 ''')

    # generate index file
    j.command(f''' bcftools index -f {j.ofile['bcf']} ''')
    return j


if __name__ == '__main__':
    backend = hb.ServiceBackend(billing_project='diverse-pop-seq-ref',
                                remote_tmpdir='gs://bge/glimpse2')  # set up backend

    # do this for each chromosome - one batch per chr
    for n in range(1, 23):

        b = hb.Batch(backend=backend,
                     name=f'GLIMPSE2-chr{n}-20KEMRI-using_GLs_and_ne')  # define batch for each chromosome

        # read in genotype likelihoods - outputs obtained from running glimpse1_gl.py
        vcf_path = f'gs://bge/glimpse1/20KEMRI/genotype_likelihoods/merged.chr{n}.vcf.gz'
        gl_vcf = b.read_input_group(**{'vcf.gz': vcf_path,
                                       'vcf.gz.csi': f'{vcf_path}.csi'})

        # a reference panel of haplotypes (only bi-allelic snps) -  files obtained from running ref_qc.py
        bcf_path = f'gs://bge/reference_panel/hgdp.tgp.trio/hgdp.tgp.chr{n}.onlysnps.bcf'
        ref_bcf = b.read_input_group(**{'bcf': bcf_path,
                                        'bcf.csi': f'{bcf_path}.csi'})

        # a fine-scale genetic map - b38
        genetic_map = b.read_input(f'gs://hgdp-1kg/phasing/maps/b38/chr{n}.b38.gmap.gz')

        # file containing the imputation chunks per chr - obtained from GLIMPSE developers (4cM)
        chunk_file = pd.read_csv(f'gs://hgdp-1kg/phasing/chunks/b38/4cM/chunks_chr{n}.txt', sep='\t', header=None)
        regions = [(irg, org) for irg, org in
                   zip(chunk_file[2], chunk_file[3])]  # grab buffered and imputation regions from txt file

        # get file sizes to allocate resources for batch jobs accordingly
        vcf_size = round(get_file_size(vcf_path))
        bcf_size = round(get_file_size(bcf_path))
        storage_impute = round(vcf_size + bcf_size + 2)
        storage_ligate = round(vcf_size * 2.5)

        results = []  # empty list for imputed chunks - STEP 2 glimpse2_phase outputs

        # for each buffered and imputation regions, do the following
        for i in range(len(regions)):

            regs = regions[i][0].split(':')[1].split('-')[0]  # start region
            rege = regions[i][0].split(':')[1].split('-')[1]  # end region

            # RUN STEP 1 - create binary reference panel
            run_binary_ref = create_binary_ref(ref_bcf, genetic_map, f'chr{n}', regs, rege, regions[i][0], regions[i][1])

            # RUN STEP 2  - Running GLIMPSE2 with genotype likelihoods that were used for GLIMPSE1
            run_glimpse2_phase = glimpse2_phase(gl_vcf, run_binary_ref.ofile, f'chr{n}', regs, rege, storage_impute)

            # append outputs to the list 'results'
            results.append(run_glimpse2_phase.ofile)

        # RUN STEP 3 - ligate chunks of the the same chromosome
        run_ligate_chunks = ligate_chunks(results, f'chr{n}', storage_ligate)

        # write out merged GLs for each chr
        b.write_output(run_ligate_chunks.ofile, f'gs://bge/glimpse2/20KEMRI/glimpse2_imputation_with_gl_and_ne/gl.ne.imputed.ligated.chr{n}')

        b.run(wait=False)  # run batch

    backend.close()  # close a batch backend







