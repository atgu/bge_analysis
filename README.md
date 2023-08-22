# bge_analysis

## Overview 

ref_qc.py - Using bcftools Filter a reference panel to only bi-allelic SNPS
glimpse1_gl.py - Using bcftools Extract the variable positions/sites from the subsetted ref panel (output from ref_qc.py); Compute genotype likelihoods (GLs) using the 
individual cram 
files; Merge the GLs for all target samples per chromosome
glimpse1_imp.py - Using GLIMPSE1 Impute and phase chromosome chunks and ligate them  back together per chromosome
glimpse2.py - Using  GLIMPSE2 Convert the reference panel into GLIMPSE2’s binary file format; Impute and phase a whole chromosome; Ligate imputed chunks of the the same 
chromosome
