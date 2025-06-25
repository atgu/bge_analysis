__author__ = 'Mary T. Yohannes and Toni Boltz'

import hailtop.fs as hfs
import hailtop.batch as hb
import subprocess
import os
import pandas as pd
import hail as hl


def glimpse2_phase(b, cram_list, ref_binary, ref_genome, chr, regs, rege, mem_type):
    j = b.new_job(name=f'run GLIMPSE2_phase - {chr}:{regs}-{rege}')  # define job

    j._preemptible = True
    j.cpu(4).memory(mem_type)
    j.image(
        'us.gcr.io/broad-dsde-methods/ckachulis/glimpse_for_wdl_pipeline:checkpointing_and_extract_num_sites')  # Docker image with GLIMPSE2 checkpointing

    # output handler
    j.declare_resource_group(ofile={
        'bcf': '{root}.bcf',
        'bcf.csi': '{root}.bcf.csi'})
  
    #mount Terra buckets to VM
    j.cloudfuse('fc-8bcb1af7-a219-4671-90f7-e5653e87c3ad', '/home/tboltz/mount_folder1/')
    j.cloudfuse('fc-d6e8ddd0-6b00-4f7c-8051-351ab37f79c1', '/home/tboltz/mount_folder2/')
    j.cloudfuse('fc-f1d2d508-6593-48ee-995d-acb46f5ec3c1', '/home/tboltz/mount_folder3/')
    j.cloudfuse('fc-9581c68b-4ad7-41d8-9b6a-cc7f6b624a62', '/home/tboltz/mount_folder4/')
    j.cloudfuse('fc-5058e99a-3e69-47a7-bc79-73aae6b19e95', '/home/tboltz/mount_folder5/')
    
    # run GLIMPSE2_phase
    j.command(f''' GLIMPSE2_phase --bam-list {cram_list} \
    --fasta {ref_genome.fasta} \
    --reference {ref_binary} \
    --threads 4 \
    --output {j.ofile['bcf']} ''')

    # generate index file
    j.command(f''' bcftools index -f {j.ofile['bcf']} ''')
    return j


def ligate_chunks(imputed_chunk_list, chr):
    j = b.new_job(name=f'ligate chunks - {chr}')  # define job
    j.cpu(4)
    j._preemptible = True
    j.image(
        'us.gcr.io/broad-dsde-methods/ckachulis/glimpse_for_wdl_pipeline:checkpointing_and_extract_num_sites') 
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
    backend = hb.ServiceBackend(billing_project='neale-pumas-bge',
                                remote_tmpdir='gs://neurogap-bge-imputed-regional',regions=['us-central1'])  # set up backend

    #this file has a list of file paths, one for each batch
    #a batch is 200 cram paths to be imputed together
    sample_batch_list = hfs.open('gs://neurogap-bge-imputed-regional/sample_paths_files/neurogap_batch_list_wpaths.txt')

    for batch in sample_batch_list.readlines():

        batch=batch.strip("\n")
        batch_name=os.path.basename(batch)

        b = hb.Batch(backend=backend, name=f'GLIMPSE2-200NeuroGAP-{batch_name}')

        cram_paths = b.read_input(batch)

        ref_genome = b.read_input_group(
            fasta='gs://gcp-public-data--broad-references/hg38/v0/Homo_sapiens_assembly38.fasta',
            fai='gs://gcp-public-data--broad-references/hg38/v0/Homo_sapiens_assembly38.fasta.fai')
        
        for n in range(1, 23):
            
            mem_type="standard" #if any CHR fails due to "Out of Memory", change this to "highmem"

            # a fine-scale genetic map - b38                                                                                               
            genetic_map = b.read_input(f'gs://hgdp-1kg/phasing/maps/b38/chr{n}.b38.gmap.gz')

            # file containing the imputation chunks per chr - obtained from GLIMPSE developers (4cM)                 
            chunk_file = pd.read_csv(f'gs://hgdp-1kg/phasing/chunks/b38/4cM/chunks_chr{n}.txt', sep='\t', header=None)
            regions = [(irg, org) for irg, org in
                       zip(chunk_file[2], chunk_file[3])]  # grab buffered and imputation regions from txt file                
        
        
            results = []  # empty list for imputed chunks
            # for each buffered and imputation regions, do the following
            for i in range(len(regions)):
            
                regs = regions[i][0].split(':')[1].split('-')[0]  # start region
                rege = regions[i][0].split(':')[1].split('-')[1]  # end region

                #get binary ref created from other script
                binary_ref = b.read_input(f'gs://neurogap-bge-imputed-regional/binary_reference/ref_chunk_chr{n}_{i}.bin')
            
                run_glimpse2_phase = glimpse2_phase(b, cram_paths, binary_ref, ref_genome, f'chr{n}', regs, rege, mem_type)

                results.append(run_glimpse2_phase.ofile)

            run_ligate_chunks = ligate_chunks(results, f'chr{n}')

            b.write_output(run_ligate_chunks.ofile, f'gs://neurogap-bge-imputed-regional/glimpse2/{batch_name}.imputed.chr{n}')

            b.run(wait=False)  # run batch

    backend.close()  # close a batch backend
