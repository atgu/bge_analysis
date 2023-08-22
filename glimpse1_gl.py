__author__ = 'Mary T. Yohannes'

### GLIMPSE1 SCRIPT #1 ###

import hailtop.batch as hb
import subprocess
import os

# In this script, you'll find how to (using bcftools):
# # Extract the variable positions/sites from the subsetted ref panel (output from ref_qc.py)
# # Compute genotype likelihoods (GLs) using the individual cram files
# # Merge the GLs for all target samples per chromosome
### Based on the tutorial probided here https://odelaneau.github.io/GLIMPSE/glimpse1/tutorial_b38.html#run_preliminaries

# Step 1) Extracting variable positions (sites) in the reference panel
def extract_sites(ref_snps):
    j = b.new_job(name="extract variable sites")  # define job
    j.image('dockerbiotools/bcftools') # Docker image with bftools and htslib
    j.cpu(4)

    # output handler
    j.declare_resource_group(ofile={
        'vcf.gz': '{root}.vcf.gz',
        'vcf.gz.csi': '{root}.vcf.gz.csi',
        'tsv.gz': '{root}.tsv.gz',
        'tsv.gz.tbi': '{root}.tsv.gz.tbi'})

    # extract sites in the reference panel
    j.command(f'''bcftools view -G -m 2 -M 2 -v snps {ref_snps['bcf']} -Oz -o {j.ofile['vcf.gz']}''')

    # generate index file
    j.command(f''' bcftools index -f {j.ofile['vcf.gz']} ''')

    # convert to tsv
    j.command(f'''bcftools query -f'%CHROM\t%POS\t%REF,%ALT\n' {j.ofile['vcf.gz']} | bgzip -c > {j.ofile['tsv.gz']}''')

    # generate index file using tabix - requires htslib
    j.command(f''' tabix -s1 -b2 -e2 {j.ofile['tsv.gz']} ''')
    return j


# Step 2) Computing genotype likelihoods (GLs) for a single individual at specific positions (NOT USING THE CRAM FILES UNTIL NOW)
def compute_gl(sample_cram, ref_sites, ref_genome, label, chr):
    j = b.new_job(name=f'compute GLs:{label}-{chr}')  # define job
    j.cpu(4)
    j.storage('10Gi')  # increase storage
    j.image('docker.io/lindonkambule/shapeit5_2023-05-05_d6ce1e2:v5.1.1')  # Docker image with bftools and htslib

    # output handler
    j.declare_resource_group(ofile={
        'vcf.gz': '{root}.vcf.gz',
        'vcf.gz.csi': '{root}.vcf.gz.csi'})

    # compute GLs
    j.command(f''' bcftools mpileup -f {ref_genome.fasta} -I -E -a 'FORMAT/DP' -T {ref_sites['vcf.gz']} -r {chr} {sample_cram.cram} -Ou | bcftools call -Aim -C alleles -T {ref_sites['tsv.gz']} -Oz -o {j.ofile['vcf.gz']}''')

    # generate index file
    j.command(f''' bcftools index -f {j.ofile['vcf.gz']} ''')
    return j


# Step 3) Merging genotype likelihoods (GLs) of multiple individuals per chr
def merge_gl(gl_list, chr):
    j = b.new_job(name=f'merge GLs-{chr}')  # define job
    j.cpu(4)
    j.image('docker.io/lindonkambule/shapeit5_2023-05-05_d6ce1e2:v5.1.1')  # Docker image with bftools and htslib

    # output handler
    j.declare_resource_group(ofile={
        'vcf.gz': '{root}.vcf.gz',
        'vcf.gz.csi': '{root}.vcf.gz.csi'})

    # the input file is a list of GLs of multiple samples (vcf.gz,vcf.gz.csi)
    # here we are getting the names/paths of those vcf.gz files
    gl_file_names = '\n'.join([f'{s["vcf.gz"]}' for s in gl_list])

    # create a text file containing the full list of vcf.gz files, one sample file per line
    j.command(f'echo "{gl_file_names}" > list_merge.txt')

    # generate a single file containing genotype likelihoods for all target individuals for a particular chr
    j.command(f'''bcftools merge -m none -r {chr} -Oz -o {j.ofile['vcf.gz']} -l list_merge.txt''')

    # generate index file
    j.command(f''' bcftools index -f {j.ofile['vcf.gz']} ''')
    return j


if __name__ == '__main__':
    backend = hb.ServiceBackend(billing_project='diverse-pop-seq-ref', remote_tmpdir='gs://bge/glimpse1') # set up backend

    # do this for each chromosome - one batch per chr
    for n in range(1, 23):

        b = hb.Batch(backend=backend, name=f'GLIMPSE1-GLs-chr{n}-20KEMRI') # define batch

        # read in a reference genome in order to align the target sample cram file - GRCH38 in our case
        ref_genome = b.read_input_group(
            fasta='gs://gcp-public-data--broad-references/hg38/v0/Homo_sapiens_assembly38.fasta',
            fai='gs://gcp-public-data--broad-references/hg38/v0/Homo_sapiens_assembly38.fasta.fai')


        # read in ref panel with only biallelic snps from the cloud - output from fer_qc.py"
        ref_onlysnps_path = f'gs://bge/reference_panel/hgdp.tgp.trio/hgdp.tgp.chr{n}.onlysnps.bcf'
        ref_onlysnps = b.read_input_group(**{'bcf': ref_onlysnps_path,
                                             'bcf.csi': f'{ref_onlysnps_path}.csi'})

        # RUN STEP 1 - extracting sites from reference panel to calculate GLs
        run_extract_sites = extract_sites(ref_onlysnps)

        # get cram file google cloud paths for target samples
        root_path = 'gs://fc-d6e8ddd0-6b00-4f7c-8051-351ab37f79c1/NeuroGAP-Psychosis_KEMRI_Psychosis_BGE/RP-1889.PDO-29432/Exome' # common root path
        sub_path = '**/**/*.cram'
        cram_paths = subprocess.check_output(['gsutil', 'ls', f'{root_path}/{sub_path}']).decode().strip('\n').split('\n') # returns a list of individual cram file paths

        # USE ONLY 20 SAMPLES FOR NOW
        cram_paths = cram_paths[0:20]

        results = []  # empty list for outputs
        for path in cram_paths:

            # only file name - without path (removed by os.path.basename) and .cram (removed by os.path.splitext) - for keeping track
            label = os.path.splitext(os.path.basename(path))[0]  # ex output: 15001708095242PK

            # read in both the .cram and .crai files
            sample_cram = b.read_input_group(
                cram=path,
                crai=f'{path}.crai')

            # RUN STEP 2 - computing GLs for a single individual at the specific sites
            run_compute_gl = compute_gl(sample_cram, run_extract_sites.ofile, ref_genome, label, f'chr{n}')

            results.append(run_compute_gl.ofile)  # append outputs to the list 'results'

        # RUN STEP 3 - merging genotype likelihoods of multiple individuals (per chromosome)
        run_merge_gl = merge_gl(results, f'chr{n}')

        # write out merged GLs for each chr
        b.write_output(run_merge_gl.ofile, f'gs://bge/glimpse1/20KEMRI/genotype_likelihoods/merged.chr{n}')

        b.run(wait=False) # run batch

    backend.close()   # close a batch backend







