__authors__ = 'Toni Boltz and Lindo Nkambule'

import hailtop.batch as hb
import subprocess
import os

def convert_to_vcf(b, batch, chrom):
    
    j = b.new_job(name=f'convert-bcf-to-vcf-{chrom}')    
    j.cpu(4)
    j.storage('10Gi')
    j.image('us.gcr.io/broad-dsde-methods/glimpse:palantir-workflows_20c9de0')

    batch_name = "neurogap_gsa_batch" + batch + ".imputed.chr" + str(chrom)
    path = "gs://neurogap-bge-imputed-regional/glimpse2/gsa/"
    full_path = path + batch_name

    bcf_batch=b.read_input_group(
        bcf=f'{full_path}.bcf',
        csi=f'{full_path}.bcf.csi')
    
    j.command(f'''bcftools view {bcf_batch['bcf']} -Ov -o {j.ofile}''')

    return j


def convertVCF2table(vcf_batch, chrom):

    j = b.new_job(name=f'convert-vcfs2table-{chrom}')

    j.cpu(4)
    j.storage('10Gi')
    j.image('us.gcr.io/broad-gatk/gatk:4.3.0.0')

    # output handler                                                                                                                 
    j.declare_resource_group(ofile={
        'tsv.gz': '{root}.tsv.gz'})

    j.command(f'''gatk VariantsToTable -V {vcf_batch} -O temp.tsv -F CHROM -F POS -F REF -F ALT -F AF -F INFO ''')
    j.command(f'''bgzip -c temp.tsv > {j.ofile['tsv.gz']}''')

    return j


def af_info_python(file_list, num_sample_list, chrom):

    p = b.new_job(name=f'fix-af-info-annot-chr{chrom}')  # define job  
    p.cpu(4)
    p.storage('10Gi')
    
    p.image('docker.io/tboltz/python-fix_af_info-script')

    i = 0
    for file in file_list:
        p.command(
            f"""
            cp {file.ofile['tsv.gz']} file_{i}.tsv.gz
        """)
        i += 1
        
    #python code below copied from: https://github.com/broadinstitute/palantir-workflows/blob/main/GlimpseImputationPipeline/Glimpse2MergeBatches.wdl#L158
    p.command(f'''
    python -c '

import functools
import pandas as pd
import glob
import os

path = "*.tsv.gz"

input_filenames = sorted(glob.glob(f"{{os.getcwd()}}/*.tsv.gz"), key=len)
    
num_samples = {num_sample_list}
    
num_batches = len(input_filenames)

def calculate_af(row):
    return sum([row[f"AF_{{i}}"] * num_samples[i] for i in range(num_batches)]) / sum(num_samples)

def calculate_info(row):
    aggregated_af = row["AF"]
    return 1 if (aggregated_af == 0.0) or (aggregated_af == 1.0) else \
        1 - (sum([(1 - row[f"INFO_{{i}}"]) * 2 * num_samples[i] * row[f"AF_{{i}}"] * (1 - row[f"AF_{{i}}"]) for i in range(num_batches)])) / \
            (2 * sum(num_samples) * aggregated_af * (1 - aggregated_af))

annotation_dfs = [pd.read_csv(input_filename, sep="\t").rename(columns={{"AF": f"AF_{{i}}", "INFO": f"INFO_{{i}}"}}) for i, input_filename in enumerate(input_filenames)]

annotations_merged = functools.reduce(lambda left, right: pd.merge(left, right, on=["CHROM", "POS", "REF", "ALT"], how="inner", validate="one_to_one"), annotation_dfs)

annotations_merged["AF"] = annotations_merged.apply(lambda row: calculate_af(row), axis=1)
annotations_merged["INFO"] = annotations_merged.apply(lambda row: calculate_info(row), axis=1)

annotations_merged.to_csv("temp_file.csv", sep="\t", columns=["CHROM", "POS", "REF", "ALT", "AF", "INFO"], header=False, index=False)
'
    ''')

    p.command(f"""
        cp temp_file.csv {p['ofile']}
    """)

    return p


def merge_batches(bcf_batch_list, annot, chrom):
    m = b.new_job(name=f'merge-batches-reannotate-{chrom}')  # define job     

    m.cpu(4)
    m.storage('10Gi')
    m.image('us.gcr.io/broad-dsde-methods/glimpse:palantir-workflows_20c9de0')

    m.declare_resource_group(ofile={
        'vcf.gz': '{root}.vcf.gz',
        'vcf.gz.csi': '{root}.vcf.gz.csi'})

    bcf_batch_names = ' '.join([f"{v['bcf']}" for v in bcf_batch_list])

    m.command(f'''bcftools merge {bcf_batch_names} -Oz -o merged_batches.vcf.gz''')

    m.command(f'''bgzip -c {annot} > {annot}.gz''')
    m.command(f'''tabix -s1 -b2 -e2  {annot}.gz''')
    
    m.command(f'''bcftools annotate -a {annot}.gz -c CHROM,POS,REF,ALT,AF,INFO -O z -o {m.ofile['vcf.gz']} merged_batches.vcf.gz ''')
    m.command(f'''bcftools index {m.ofile['vcf.gz']}''')

    return m
    

if __name__ == '__main__':

    backend = hb.ServiceBackend(billing_project='neale-pumas-bge',
                                remote_tmpdir='gs://neurogap-bge-imputed-regional',regions=['us-central1'])  # set up backend                                                            
    b = hb.Batch(backend=backend, name=f'merge-bcfs') # define batch   
    
    sample_batches = ['00', '01', '02', '03', '04']

    num_samples = [200, 200, 200, 200, 54]
    
    for n in range(1, 22):
        
        run_convert_to_vcf = [convert_to_vcf(b, batch,  f'{n}') for batch in sample_batches] 
        run_vcf2table = [convertVCF2table(i.ofile,f'chr{n}') for i in run_convert_to_vcf]

        # run python script to create aggregated annotation files per chromosome across all batches
        annot_file = af_info_python(run_vcf2table, num_samples, f'{n}')

        batch_list=[b.read_input_group(
                bcf=f'gs://neurogap-bge-imputed-regional/glimpse2/gsa/neurogap_gsa_batch{i}.imputed.chr{n}.bcf',                                                
                csi=f'gs://neurogap-bge-imputed-regional/glimpse2/gsa/neurogap_gsa_batch{i}.imputed.chr{n}.bcf.csi') for i in sample_batches]
        
        run_merge_batches = merge_batches(batch_list, annot_file.ofile, f'chr{n}')

        b.write_output(run_merge_batches.ofile, f'gs://neurogap-bge-imputed-regional/glimpse2/merged/gsa_merged_chr{n}')

    b.run(wait=False) # run batch

    backend.close()   # close a batch backend
