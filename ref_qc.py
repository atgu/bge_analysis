__author__ = 'Mary T. Yohannes'

### SCRIPT #1 ###

# In this script, you'll find how to (using bcftools):
# # Filter a reference panel to only bi-allelic SNPS
### Based on the tutorial probided here https://odelaneau.github.io/GLIMPSE/glimpse1/tutorial_b38.html#run_preliminaries

import hailtop.batch as hb
import subprocess
import os

# perform basic QC on the reference panel - only keep biallelic snps
def only_snps(ref_panel, chr):
    j = b.new_job(name=f'{chr}-biallelic snps only')  # define job
    j.cpu(4)
    j.image('docker.io/lindonkambule/shapeit5_2023-05-05_d6ce1e2:v5.1.1') # Docker image with bftools

    # output handler
    j.declare_resource_group(ofile={
        'bcf': '{root}.bcf',
        'bcf.csi': '{root}.bcf.csi'})

    # only keep SNPs - remove multiallelic records since GLIMPSE only handles bi-allelic SNPs and bcftools does not perform a good job at calling indels
    j.command(f'''bcftools norm -m -any {ref_panel['bcf']} -Ou --threads 4 | bcftools view -m 2 -M 2 -v snps --threads 4 -Ob -o {j.ofile['bcf']} ''')

    # generate index file
    j.command(f''' bcftools index -f {j.ofile['bcf']} --threads 4 ''')
    return j

if __name__ == '__main__':
    backend = hb.ServiceBackend(billing_project='diverse-pop-seq-ref', remote_tmpdir='gs://bge') # set up backend

    b = hb.Batch(backend=backend, name=f'subset ref to biallelic snps only')  # define batch

    # do this for each chromosome - one batch per chr
    for n in range(1, 23):

        # read in a reference panel based on chromosome
        ref_path = f'gs://gcp-public-data--gnomad/resources/hgdp_1kg/phased_haplotypes_v2/hgdp1kgp_chr{n}.filtered.SNV_INDEL.phased.shapeit5.bcf'
        ref_panel = b.read_input_group(**{'bcf': ref_path,
                                          'bcf.csi': f'{ref_path}.csi'})

        # run basic QC on the reference panel so we are keeping only keep bi-allelic snps
        run_only_snps = only_snps(ref_panel, f'chr{n}')

        # write out the ref chromosome with only bi-allelic snps to a google cloud bucket - to use in both GLIMPSE1/2
        b.write_output(run_only_snps.ofile, f'gs://bge/reference_panel/hgdp.tgp.trio/hgdp.tgp.chr{n}.onlysnps')

    b.run(wait=False)  # run batch

    backend.close()  # close a batch backend



