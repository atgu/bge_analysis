__author__ = 'Mary T. Yohannes'

### SCRIPT #pre-1 ###

# In this script, you'll find how to (using bcftools):
# # produce a reference panel with multiple minor allele count (MAC) >= 1,2,3,4,5 cutoffs

import hailtop.batch as hb
import subprocess
import os


def filter_mac(input_bcf, chr, i):
    j = b.new_job(name=f'{chr}')  # define job
    j.cpu(8)
    j.image('docker.io/lindonkambule/shapeit5_2023-03-23_a4a1818:latest')  # Docker image with bftools

    j.declare_resource_group(ofile={
        'bcf': '{root}.bcf',
        'bcf.csi': '{root}.bcf.csi'})

    # check number of variants before filtering out singletons (MAC<2)
    j.command(f"""
    echo Initial number of variants before QC
    bcftools query -f '%POS\n' {input_bcf['bcf']} | wc -l
""")

    j.command(f"""
    echo 1. QC
    bcftools view -i'MAC>={i}' --output {j.ofile['bcf']} {input_bcf['bcf']} 

    echo 2. Indexing
    bcftools index {j.ofile['bcf']} --output {j.ofile['bcf.csi']} --threads {8}

    echo 3. Number of variants after  QC
    bcftools query -f '%POS\n' {j.ofile['bcf']} | wc -l

    echo 3. Writing file
""")
    return j


if __name__ == '__main__':
    backend = hb.ServiceBackend(billing_project='diverse-pop-seq-ref', remote_tmpdir='gs://bge')  # set up backend

    for i in range(5, 6):
        b = hb.Batch(backend=backend, name=f'filter ref on MAC>={i}')  # define batch

        for n in range(1, 23):

            # read in a reference panel based on chromosome - this is the QC'ed ref panel that only has biallelic snps (output from ref_qc.py script)
            ref_path = f'gs://bge/reference_panel/hgdp.tgp.trio/hgdp.tgp.chr{n}.onlysnps.bcf'
            ref_panel = b.read_input_group(**{'bcf': ref_path,
                                              'bcf.csi': f'{ref_path}.csi'})

            # run MAC filtering
            run_filter = filter_mac(ref_panel, f'chr{n}', f'{i}')

            # write out the ref chromosome with only bi-allelic snps to a google cloud bucket - to use in both GLIMPSE1/2
            b.write_output(run_filter.ofile, f'gs://bge/reference_panel/hgdp.tgp.trio/mac_cutoff/mac_{i}/hgdp.tgp.onlysnps.chr{n}.mac{i}')

    b.run(wait=False)  # run batch

    backend.close()  # close a batch backend
