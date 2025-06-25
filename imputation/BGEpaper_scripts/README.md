Use `glimpse2_cram_chkpt_batches.py` to run batches of 200 individuals at a time (to minimize costs) through the GLIMPSE2 imputation software on the Broad's Hail Batch service. Each batch of 200 individuals will be written to its own .bcf file. 

Use `fix_annotations_and_merge_batches.py` to update the INFO scores and allele frequencies after re-merging the batches post-imputation. This script is a Hail Batch implementation of code described here: 
https://github.com/broadinstitute/palantir-workflows/blob/main/GlimpseImputationPipeline/README.md#glimpse2mergebatches 

