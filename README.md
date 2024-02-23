# bge_analysis

## Overview 

GitHub Repo for BGE analysis 

`glimpse2_cram_chkpt_batches.py`

Run batches of 200 individuals at a time through the GLIMPSE2 imputation software on the Broad's Hail Batch service. Each batch of 200 individuals will be written to its own .bcf file. 

`fix_annotations_and_merge_batches.py`

This script is a Hail Batch implementation of code described here: 
https://github.com/broadinstitute/palantir-workflows/blob/main/GlimpseImputationPipeline/README.md#glimpse2mergebatches 

The .bcf files outputted from the glimpse2_cram_chkpt_batches.py script are first converted to .vcf and then inputted into the GATK VariantsToTable function. These tables are then inputted into a python script which re-computes the AF and INFO scores for all batches combined. The output from the python script is an annotation file that is then used in the final step. The original .bcf files are merged with `bcftools` and then reannotated with the new annotation file. 

