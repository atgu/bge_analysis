__author__ = 'Mary T. Yohannes and Toni Boltz'

import hailtop.batch as hb
import subprocess
import os
import pandas as pd
import hail as hl


def create_binary_ref(ref_bcf, genetic_map, chr, regs, rege, irg, org):
    j = b.new_job(name=f'create binary ref panel - {chr}:{regs}-{rege}')  # define job
    j.cpu(4)
    j.storage('10Gi')
    j.image(
        'docker.io/imary116/glimpse2:with-bcftools-and-updated-info-score')  # Docker image with GLIMPSE2 with INFO score corrected

    j.command(f''' GLIMPSE2_split_reference --reference {ref_bcf['bcf']} \
        --map {genetic_map} \
        --input-region {irg} \
        --output-region {org} \
        --threads 4 \
        --output binary_reference_panel 

mv binary_reference_panel_{chr}_{regs}_{rege}.bin {j.ofile}''')
    return j

if __name__ == '__main__':
    backend = hb.ServiceBackend(billing_project='diverse-pop-seq-ref',
                                remote_tmpdir='gs://neurogap-bge-imputed')  # set up backend

    b = hb.Batch(backend=backend, name=f'create-binary-reference')

    ref_genome = b.read_input_group(
    fasta='gs://gcp-public-data--broad-references/hg38/v0/Homo_sapiens_assembly38.fasta',
    fai='gs://gcp-public-data--broad-references/hg38/v0/Homo_sapiens_assembly38.fasta.fai')

    for n in range(1, 23):

        # a reference panel of haplotypes (only bi-allelic snps)                                                    
        bcf_path = f'gs://bge/reference_panel/hgdp.tgp.trio/mac_cutoff/mac_2/hgdp.tgp.onlysnps.chr{n}.mac2.bcf'
        ref_bcf = b.read_input_group(**{'bcf': bcf_path,
                                        'bcf.csi': f'{bcf_path}.csi'})

        # a fine-scale genetic map - b38                                                                                                                                       
        genetic_map = b.read_input(f'gs://hgdp-1kg/phasing/maps/b38/chr{n}.b38.gmap.gz')

        # file containing the imputation chunks per chr - obtained from GLIMPSE developers (4cM)                                                                               
        chunk_file = pd.read_csv(f'gs://hgdp-1kg/phasing/chunks/b38/4cM/chunks_chr{n}.txt', sep='\t', header=None)
        regions = [(irg, org) for irg, org in
                   zip(chunk_file[2], chunk_file[3])]  # grab buffered and imputation regions from txt file                
        
        # for each buffered and imputation regions, do the following
        for i in range(len(regions)):

            regs = regions[i][0].split(':')[1].split('-')[0]  # start region
            rege = regions[i][0].split(':')[1].split('-')[1]  # end region

            run_binary_ref = create_binary_ref(ref_bcf, genetic_map, f'chr{n}', regs, rege, regions[i][0], regions[i][1])

            b.write_output(run_binary_ref.ofile, f'gs://neurogap-bge-imputed/binary_reference/ref_chunk_chr{n}_{i}.bin')

        b.run(wait=False)  # run batch

    backend.close()  # close a batch backend
