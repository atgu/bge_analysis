Use the scripts provided in this directory to complete the Wiki tutorial. 

Steps for PRS (example provided here using Ethiopian NeuroGAP cohort): 
1. `PRS_QC.ipynb` script to quality control BGEN files
2. `fmt_files.sh` to get all files formatted for PRS-CS
3. `run_prscs_forWiki.sh` to run the PRS-CS software (see https://github.com/getian107/PRScs)
4. `weights.sh` to concatenate all the per-chromosome SNP weights and compute scores for the cohort
5. `prs_metrics_forWiki.R` to compute predictiability metrics (including Nagelkerke R2, liability adjustment, and AUC) and produce plots

