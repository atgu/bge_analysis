__author__ = 'Mary T. Yohannes'

### GLIMPSE1 SCRIPT #2 ###

import hailtop.batch as hb
import subprocess
import os
import hail as hl
import pandas as pd

# In this script, you'll find how to (using GLIMPSE):
# # Impute and phase chromosome chunks
# # Ligate the imputed chunks back together per chromosome
### Based on the tutorial probided here https://odelaneau.github.io/GLIMPSE/glimpse1/tutorial_b38.html#run_preliminaries


# get file size to allocate resources for batch jobs accordingly
def get_file_size(file):
    file_info = hl.utils.hadoop_stat(file)
    size_bytes = file_info['size_bytes']
    size_gigs = size_bytes / (1024 * 1024 * 1024)
    return size_gigs

# Step 1) Impute and phase a whole chromosome using chunks
def impute_chunk(gl_vcf, ref_bcf, genetic_map, irg, org, storage_size):
    j = b.new_job(name=f'impute chunk-{org}')  # define job
    j.cpu(4)
    j.storage(f'{storage_size}Gi')
    j.image('docker.io/simrub/glimpse:v1.1.1-c27e90d_20210521')  # Docker image with GLIMPSE and bcftools

    # output handler
    j.declare_resource_group(ofile={
        'bcf': '{root}.bcf',
        'bcf.csi': '{root}.bcf.csi',
        'log': '{root}.log'})

    j.command(f''' GLIMPSE_phase_v1.1.1 --input {gl_vcf['vcf.gz']} \
    --reference {ref_bcf['bcf']} \
    --map {genetic_map} \
    --input-region {irg} \
    --output-region {org} \
    --output {j.ofile['bcf']} \
    --log {j.ofile['log']} ''') # --thread {ncpu - 2} \

    # generate index file
    j.command(f''' bcftools index -f {j.ofile['bcf']} ''')
    return j

# Step 2) Ligate the imputed chunks together
def ligate_chunks(imputed_chunk_list, chr, storage_size):
    j = b.new_job(name=f'ligate chunks-{chr}')  # define job
    j.cpu(4)
    j.storage(f'{storage_size}Gi')
    j.image('docker.io/simrub/glimpse:v1.1.1-c27e90d_20210521')  # Docker image with GLIMPSE and bcftools

    # output handler
    j.declare_resource_group(ofile={
        'bcf': '{root}.bcf',
        'bcf.csi': '{root}.bcf.csi',
        'log': '{root}.log'})

    # the input file is a list of imputed chunks
    # here we are getting the names/paths of those chunks
    imp_chunk_bcf_names = '\n'.join([f'{c["bcf"]}' for c in imputed_chunk_list])

    # create a text file containing the full list of bcf files, one imputed chunk file per line
    j.command(f'echo "{imp_chunk_bcf_names}" > chunk_list.txt')

    j.command(f''' GLIMPSE_ligate_v1.1.1 --input chunk_list.txt --output {j.ofile['bcf']} --log {j.ofile['log']} ''')

    # generate index file
    j.command(f''' bcftools index -f {j.ofile['bcf']} ''')
    return j

if __name__ == '__main__':
    backend = hb.ServiceBackend(billing_project='diverse-pop-seq-ref', remote_tmpdir='gs://bge/glimpse1') # set up backend

    # for each chromosome, do the following
    for n in range(1, 23):

        b = hb.Batch(backend=backend, name=f'GLIMPSE1-impute&ligate-chr{n}')  # define batch for each chromosome

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
        regions = [(irg, org) for irg, org in zip(chunk_file[2], chunk_file[3])] # grab buffered and imputation regions from txt file

        # get file sizes to allocate resources for batch jobs accordingly
        vcf_size = round(get_file_size(vcf_path))
        bcf_size = round(get_file_size(bcf_path))
        storage_impute = round(vcf_size + bcf_size + 2)
        storage_ligate = round(vcf_size * 2.5)

        # empty list for imputed chunks
        impute_results = []

        # for each buffered and imputation regions, do the following
        for i in range(len(regions)):

            # RUN STEP 1 - impute the chunk with the specified region
            run_impute_chunk = impute_chunk(gl_vcf, ref_bcf, genetic_map, regions[i][0], regions[i][1], storage_impute)

            impute_results.append(run_impute_chunk.ofile) # append imputed chunk to output list

        # RUN STEP 2 - ligate imputed chunks (per chromosome)
        run_ligate_chunks = ligate_chunks(impute_results, f'chr{n}', storage_ligate)

        # write out merged imputed chunks for each chr
        b.write_output(run_ligate_chunks.ofile, f'gs://bge/glimpse1/20KEMRI/glimpse1_imputation/chr{n}.imputed.merged')

        b.run(wait=False)  # run batch for each chr

    backend.close()  # close a batch backend