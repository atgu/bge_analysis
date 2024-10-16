#!/bin/bash

#concatenate all per-chromsome SNP weights to one file 
cat AAU_PRSCS_pst_eff_a1_b0.5_phiauto_chr* > AAU_weights.txt

#use plink to get scores of the AAU cohort
plink2 --bfile AAU_prs_qc2 --score AAU_weights.txt cols=+scoresums ignore-dup-ids 2 4 6  --out AAU_scores

