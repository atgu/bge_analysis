# GitHub Repo for BGE analysis 
---

## Overview 

BGE Analysis paper now on BioRXiv! 
https://www.biorxiv.org/content/10.1101/2024.09.06.611689v1

---

For Broad Institute researchers looking to impute BGE data using GLIMPSE2, please see the `imputation` directory, which contains the following sub-directories:

- `glimpse_hail_batch` includes the fully documented and fully automated code for running autosomes and X chromosome imputation using GLIMPSE2 on BGE data, it assumes input files are stored in a Terra or Google Cloud bucket. For anyone working within the Broad Institute and with access to Hail Batch, these scripts should be the starting point for imputation.
  
  The following two sub-directories are included as archival references of the imputation work done during manuscript development. They are not recommended for general use, as they lack the automation and robust documentation of the current `glimpse_hail_batch` implementation:  
  
- `BGEpaper_scripts` includes the scripts used for the imputation of the data in the BGE paper. These scripts perform the same tasks as the `glimpse_hail_batch` scripts but requires more onus on the user for providing input paths and pre-splitting batches into groups of 200 individuals.
- `glimpse_comparison` includes the initial scripts used in implementing GLIMPSE2 on Hail Batch and testing the GLIMPSE parameters on a subset of BGE data.

For those outside of the Broad Institute who would like to run imputation using GLIMPSE2 on BGE data, Hail Batch is unfortunately not an option that is currently available. We recommend following the WDL provided here instead: https://github.com/broadinstitute/palantir-workflows/tree/main/GlimpseImputationPipeline 

---

The `concordance` directory contains scripts used in computing non-reference concordance and aggregated R2 metrics of BGE imputed SNPs with Global Screening Array chip data used as ground truth, as shown in the manuscript. 

---

The `exome_qc` directory contains python and R scripts for generating the QC metrics of BGE data and corresponding plots, as shown in the manuscript. 

---

For calling CNVs on BGE data, please see the repository here: https://github.com/broadinstitute/gatk/tree/4.1.0.0/scripts/cnv_wdl/germline 


