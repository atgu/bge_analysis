# bge_analysis

## Overview 

- `ref_qc.py` - using bcftools:
  - filter a reference panel to only bi-allelic SNPS 

- `glimpse1_gl.py` - using bcftools: 
  - extract variable positions/sites from the subsetted ref panel (output from `ref_qc.py`)  
  - compute genotype likelihoods (GLs) using individual cram files from NeuroGAP 
  - merge the GLs for all target samples per chromosome

- `glimpse1_imp.py` - using GLIMPSE1:
  - impute and phase chromosome chunks and ligate them  back together per chromosome using GLs 

- `glimpse2.py` - using  GLIMPSE2:
  - convert the reference panel in to GLIMPSE2’s binary file format
  - impute and phase per chromosome
  - ligate imputed chunks of the the same chromosome
  
- `glimpse2_with_gl.py`:
  - similar to `glimpse2.py` but use the GLs that were generated for GLIMPSE1 as input (`--input-gl` argument instead of `--bam-list` for `GLIMPSE2_phase`) 
  - suspect that the internal calling of the GLs in GLIMPSE2 does not work as well as bcftools/GATK GLs, e.g. in high-coverage regions 
  
- `glimpse2_with_gl_and_ne.py`: 
  - similar to `glimpse2_with_gl.py` but add `--ne 20000` argument in addition to using GLs 
  - basically the GLIMPSE1 model
  

