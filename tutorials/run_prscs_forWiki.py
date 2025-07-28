#!/bin/bash

sst="gs://neurogap-bge-imputed-regional/forPRSWiki/formatted_PGC3_SCZ.txt"
N=170114

#the LD reference was downloaded from the PRS-CS github repo, see here: https://github.com/getian107/PRScs
#specifically the AFR 1KG reference panel: ldblk_1kg_afr.tar.gz
refdir="gs://neurogap-bge-imputed-regional/toni/prscs/ldblk_1kg_afr"

for chrom in range(1, 23):
do
    bfile_path = "gs://neurogap-bge-imputed-regional/toni/forPRSWiki/Uganda_prs_chr"${chrom}

    python3 PRScs.py \
    --ref_dir=${ref_dir} \
    --bim_prefix=${bfile_path} \
    --sst_file=${sst} \
    --n_gwas=${N} \
    --out_dir=tmp_prscs_output/Uganda_PRSCS \
    --chrom=${chrom}

done

