# bge_analysis

## Overview 

GitHub Repo for BGE analysis 

BGE Analysis paper now on BioRXiv: https://www.biorxiv.org/content/10.1101/2024.09.06.611689v1

The `exome-qc` directory contains python and R scripts for generating the QC metrics of BGE data and corresponding plots as shown in the manuscript. 

The `imputation` directory contains a sub-folder called `glimpse_comparison`, which are the initial scripts used in implementing GLIMPSE on Hail Batch and testing the GLIMPSE parameters on a subset of BGE data. The final scripts used in the BGE imputation pipeline include: 

- `glimpse2_cram_chkpt_batches.py`

Run batches of 200 individuals at a time (to minimize costs) through the GLIMPSE2 imputation software on the Broad's Hail Batch service. Each batch of 200 individuals will be written to its own .bcf file. 

  - `fix_annotations_and_merge_batches.py`

This script is a Hail Batch implementation of code described here: 
https://github.com/broadinstitute/palantir-workflows/blob/main/GlimpseImputationPipeline/README.md#glimpse2mergebatches 

The `concordance` directory contains scripts used in computing accuracy metrics of BGE imputed SNPs with Global Screening Array data used as ground truth, as shown in the manuscript. 

