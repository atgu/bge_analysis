__authors__ = 'Toni Boltz and Lindo Nkambule'

import hailtop.batch as hb
import hail as hl
import subprocess
import os

def get_file_size(file):
    file_info = hl.utils.hadoop_stat(file)
    size_bytes = file_info['size_bytes']
    size_gigs = size_bytes / (1024 * 1024 * 1024)
    return size_gigs

def convertBCF2table(bcf_batch, chrom, storage):

    j = b.new_job(name=f'convert-bcfs2table-{chrom}')

    j.memory(f'{storage}Gi')
    j.storage(f'{storage}Gi')
    j.image('us.gcr.io/broad-dsde-methods/glimpse:palantir-workflows_20c9de0')

    # output handler                                                                                                                 
    j.declare_resource_group(ofile={
        'tsv.gz': '{root}.tsv.gz'})

    j.command(f'''bcftools query -f '%CHROM %POS %REF %ALT %AF %INFO/INFO\n' {bcf_batch['bcf']} > temp.tsv''')
    j.command(f''' echo -e "CHROM POS REF ALT AF INFO" > temp_wheader.tsv && cat temp.tsv >> temp_wheader.tsv ''')
    j.command(f'''bgzip -c temp_wheader.tsv > {j.ofile['tsv.gz']}''')

    return j


def af_info_python(file_list, num_sample_list, chrom):

    p = b.new_job(name=f'fix-af-info-annot-chr{chrom}')  # define job  
    p.memory("highmem")
    p.storage('1Gi')
    
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

annotation_dfs = [pd.read_csv(input_filename, sep=" ").rename(columns={{"AF": f"AF_{{i}}", "INFO": f"INFO_{{i}}"}}) for i, input_filename in enumerate(input_filenames)]

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


def merge_batches(bcf_batch_list, annot, chrom, storage):
    m = b.new_job(name=f'merge-batches-reannotate-{chrom}')  # define job     

    m.memory(f'{storage}Gi')
    m.storage(f'{storage}Gi')
    m.image('us.gcr.io/broad-dsde-methods/glimpse:palantir-workflows_20c9de0')

    m.declare_resource_group(ofile={
        'bcf': '{root}.bcf',
        'bcf.csi': '{root}.bcf.csi'})

    bcf_batch_names = ' '.join([f"{v['bcf']}" for v in bcf_batch_list])

    m.command(f'''bcftools merge {bcf_batch_names} -Ob -o /io/merged_batches.bcf''')

    m.command(f'''bgzip -c {annot} > {annot}.gz''')
    m.command(f'''tabix -s1 -b2 -e2  {annot}.gz''')
    
    m.command(f'''bcftools annotate -a {annot}.gz -c CHROM,POS,REF,ALT,AF,INFO -Ob -o {m.ofile['bcf']} /io/merged_batches.bcf ''')
    m.command(f'''bcftools index {m.ofile['bcf']}''')

    return m
    

if __name__ == '__main__':

    backend = hb.ServiceBackend(billing_project='neale-pumas-bge',
                                remote_tmpdir='gs://neurogap-bge-imputed-regional',regions=['us-central1'])  # set up backend                                                            
    b = hb.Batch(backend=backend, name=f'NGAP-wave2-merge-bcfs') # define batch   
    
    sample_batches = ['00', '01', '02', '03', '04', '05', '06', '07', '08', '09', '10', '11', '12', '13', '14', '15', '16', '17']

    num_samples = [200, 200, 200, 200, 200, 200, 200, 200, 200, 200, 200, 200, 200, 200, 200, 200, 200, 202]
    
    for n in range(1, 23):

        batch_list=[b.read_input_group(
            bcf=f'gs://neurogap-bge-imputed-regional/glimpse2/neurogap_wave2_batch{i}.imputed.chr{n}.bcf',
            csi=f'gs://neurogap-bge-imputed-regional/glimpse2/neurogap_wave2_batch{i}.imputed.chr{n}.bcf.csi') for i in sample_batches]

        batch_list_names = [f'gs://neurogap-bge-imputed-regional/glimpse2/neurogap_wave2_batch{i}.imputed.chr{n}.bcf' for i in sample_batches]
        
        bcf_sizes = [round(get_file_size(batch))for batch in batch_list_names]
        storages = [round(size) for size in bcf_sizes]

        run_bcf2table = [convertBCF2table(batch, f'chr{n}', storage) for batch, storage in zip(batch_list, storages)]
        
        # run python script to create aggregated annotation files per chromosome across all batches
        annot_file = af_info_python(run_bcf2table, num_samples, f'{n}')

        total_storage = (sum(storages))*10
        
        run_merge_batches = merge_batches(batch_list, annot_file.ofile, f'chr{n}', total_storage)

        b.write_output(run_merge_batches.ofile, f'gs://neurogap-bge-imputed-regional/glimpse2/merged_bcfs/neurogap_wave2_merged_chr{n}')

    b.run(wait=False) # run batch

    backend.close()   # close a batch backend
