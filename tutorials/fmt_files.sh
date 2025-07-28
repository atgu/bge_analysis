#!/bin/bash

#plink commands to prepare genetics files for PRS-CS

#remove genotypes with >5% missingness
plink2 --bfile Uganda --geno 0.05 --make-bed --out Uganda_prs_qc1

#remove samples with >5% missingness
plink2 --bfile Uganda_prs_qc1 --mind 0.05 --make-bed --out Uganda_prs_qc2

#get per-chromosome plink files for PRS-CS
for c in {1..22};
do

    ~/Desktop/apps/plink2 --bfile Uganda_prs_qc2 --chr ${c} --make-bed --out Uganda_prs_chr${c}

done


###########
#downloading and formatting GWAS summary statistics 

#from PGC Downloads page (https://pgc.unc.edu/for-researchers/download-results/)
#downloaded this file:    https://figshare.com/articles/dataset/scz2022/19426775?file=34517861
mv ~/Downloads/PGC3_SCZ_wave3.primary.autosome.public.v3.vcf.tsv.gz ./


gunzip PGC3_SCZ_wave3.primary.autosome.public.v3.vcf.tsv.gz

#remove headers from original .tsv file which are lines that begin with "##"
grep -v '^##' PGC3_SCZ_wave3.primary.autosome.public.v3.vcf.tsv > PGC3_SCZ_wave3.primary.autosome.public.v3.vcf_noheader.tsv

#get the median N
cut -f16 PGC3_SCZ_wave3.primary.autosome.public.v3.vcf_noheader.tsv | grep -v 'NEFFDIV2' | sort -n | awk '{a[NR]=$1} END {if (NR%2==1) {print a[(NR+1)/2]} else {print (a[NR/2] + a[NR/2+1])/2}}'
#comes out to 85057, times by 2 to get Neff = 170114

#format for PRS-CS according to https://github.com/getian107/PRScs
cut -f2,4,5,9,10 PGC3_SCZ_wave3.primary.autosome.public.v3.vcf_noheader.tsv > formatted_PGC3_SCZ_init.txt
#make sure every row has 5 fields
awk 'NR == 1 || NF == 5'  formatted_PGC3_SCZ_init.txt > formatted_PGC3_SCZ.txt

