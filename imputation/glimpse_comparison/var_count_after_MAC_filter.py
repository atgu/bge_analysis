__author__ = 'Mary T. Yohannes'

### SCRIPT #pre-1 ###

# In this script, you'll find how to (using bcftools):
# # count how many variants there are in the reference (HGDP+1kGP) that have minor allele count (MAC) >= 1,2,3,4,5
import hailtop.batch as hb
import subprocess
import os


def filter_mac(input_bcf, chr, i):
    j = b.new_job(name=f'{chr}')  # define job
    j.cpu(8)
    j.image('docker.io/lindonkambule/shapeit5_2023-03-23_a4a1818:latest')  # Docker image with bftools

    j.declare_resource_group(ofile={
        'bcf': '{root}.bcf',
        'bcf.csi': '{root}.bcf.csi',
        'txt': '{root}.txt'})

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
    
    echo 3. Writing to file
    
    echo {i} {chr} $(bcftools query -f '%POS\n' {input_bcf['bcf']} | wc -l) $(bcftools query -f '%POS\n' {j.ofile['bcf']} | wc -l) > {j.ofile['txt']} 

""")
    return j


def merge(results):
    j = b.new_job(name='merge_results')
    j.image('docker.io/hailgenetics/hail:0.2.124-py3.11')  # docker image with Hail and Python
    j.storage('20Gi')
    j.cpu(4)

    if results:
        delimiter = "', '"
        j.command(f'''
python3 -c "
import hail as hl
ht = hl.import_table(['{delimiter.join(results)}'], impute=True, no_header=True)
ht.export('{j.ofile}')"
''')
    return j


if __name__ == '__main__':
    backend = hb.ServiceBackend(billing_project='diverse-pop-seq-ref', remote_tmpdir='gs://bge')  # set up backend

    for i in range(1, 6):
        b = hb.Batch(backend=backend, name=f'count variant MAC>={i}')  # define batch

        # do this for each chromosome - one job per chr
        r = []  # for results

        for n in range(1, 23):
            # read in a reference panel based on chromosome
            ref_path = f'gs://gcp-public-data--gnomad/resources/hgdp_1kg/phased_haplotypes_v2/hgdp1kgp_chr{n}.filtered.SNV_INDEL.phased.shapeit5.bcf'
            ref_panel = b.read_input_group(**{'bcf': ref_path,
                                              'bcf.csi': f'{ref_path}.csi'})

            # run MAC filtering and count variants before and after
            run_filter = filter_mac(ref_panel, f'chr{n}', f'{i}')

            r.append(run_filter.ofile.txt)

        run_merge = merge(r)

        # write out the ref chromosome with only bi-allelic snps to a google cloud bucket - to use in both GLIMPSE1/2
        b.write_output(run_merge.ofile, f'gs://bge/reference_panel/hgdp.tgp.trio/hgdp.tgp.mac{i}.txt')

        b.run(wait=False)  # run batch

    backend.close()  # close a batch backend
