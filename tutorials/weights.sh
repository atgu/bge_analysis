#!/bin/bash

#concatenate all per-chromsome SNP weights to one file 
cat Uganda_PRSCS_pst_eff_a1_b0.5_phiauto_chr* > Uganda_weights.txt

#use plink to get scores of the AAU cohort
plink2 --bfile Uganda_prs_qc2 --score Uganda_weights.txt cols=+scoresums ignore-dup-ids 2 4 6  --out Uganda_scores

