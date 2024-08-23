
# run sample qc
# https://broadinstitute.github.io/gnomad_methods/api_reference/sample_qc/filtering.html#gnomad.sample_qc.filtering.compute_stratified_metrics_filter
import hail as hl
import gnomad
from gnomad.sample_qc.filtering import compute_stratified_metrics_filter
hl.init(driver_cores=8, worker_memory='highmem')

MT = 'bge_wave1_variant_qc.mt' # output from step 1
SAMPLE_QC_OUT = 'pumas_bge_wave1_sample_qc.tsv.bgz'
SAMPLE_QC_OUT_FAILED = 'pumas_bge_wave1_samples_failed_prefilters.tsv.bgz'
PASSING_SAMPLES = 'pumas_bge_wave1_sample_qc_passing_samples.tsv'
PASS_FAIL_OUT = 'pumas_bge_wave1_sample_qc_pass_fail_metrics.tsv'
# load pumas
mt = hl.read_matrix_table(MT)

# filter to pumas samples
meta_bge = hl.import_table('bge_wave1_manuscript_samples.tsv', key="SAMPLE_ALIAS")
mt = mt.filter_cols(hl.is_defined(meta_bge[mt.col_key]))

# remove high chim/contam samples
samples = hl.import_table('bge_manuscript_filtered_samples_chim_contam_5_exome_90_wgs_1x.tsv', key='s')
mt = mt.filter_cols(hl.is_defined(samples[mt.col_key]))

# filter to passing ancestry assignment
ancestry = hl.import_table('pumas_bge_samples_pcs_rf_prob.txt', key='s', impute=True)
ancestry = ancestry.filter(ancestry.pass_70_per_threshold=="pass")
mt = mt.filter_cols(hl.is_defined(ancestry[mt.col_key]))

# create cohort x ancestry flag
ancestry = ancestry.annotate(cohort_pop = ancestry.TERRA_WORKSPACE + "_" + ancestry.Population)

# read in variant qc'ed mt
sample_qc_ht = hl.sample_qc(mt).cols().flatten().key_by('s')
ht = sample_qc_ht.annotate(qc_pop = ancestry[sample_qc_ht.s].cohort_pop)

# save all sample qc
sample_qc_ht.export(output=SAMPLE_QC_OUT)

filtering_qc_metrics = ['sample_qc.r_ti_tv','sample_qc.n_singleton','sample_qc.n_insertion','sample_qc.n_deletion','sample_qc.n_transition','sample_qc.n_transversion','sample_qc.r_het_hom_var','sample_qc.r_insertion_deletion']
stratified_metrics_ht = compute_stratified_metrics_filter(
                ht,
                qc_metrics={metric: ht[metric] for metric in filtering_qc_metrics},
                strata={"qc_pop": ht.qc_pop},
            )
passing_sample_qc = stratified_metrics_ht.filter((stratified_metrics_ht['fail_sample_qc.r_ti_tv']==False) & 
                                        (stratified_metrics_ht['fail_sample_qc.n_singleton']==False) & 
                                        (stratified_metrics_ht['fail_sample_qc.n_insertion']==False) & 
                                        (stratified_metrics_ht['fail_sample_qc.n_deletion']==False) & 
                                        (stratified_metrics_ht['fail_sample_qc.n_transition']==False) & 
                                        (stratified_metrics_ht['fail_sample_qc.n_transversion']==False) & 
                                        (stratified_metrics_ht['fail_sample_qc.r_het_hom_var']==False) & 
                                        (stratified_metrics_ht['fail_sample_qc.r_insertion_deletion']==False)
                                            )
mt = mt.filter_cols(hl.is_defined(passing_sample_qc[mt.col_key]))
passing = passing_sample_qc.select_globals()
passing.export(PASSING_SAMPLES)

stratified_metrics_ht.select_globals().export(PASS_FAIL_OUT)

#####
## calculate sample qc metrics for samples that failed prefilters or ancestry assingment

mt = hl.read_matrix_table(MT)

# filter to pumas samples
meta_bge = hl.import_table('bge_wave1_manuscript_samples.tsv', key="SAMPLE_ALIAS")
mt = mt.filter_cols(hl.is_defined(meta_bge[mt.col_key]))

# remove samples with qc already calculated on them
done = hl.import_table(SAMPLE_QC_OUT, key='s')
mt = mt.filter_cols(~hl.is_defined(done[mt.col_key]))

# save all sample qc
sample_qc_ht = hl.sample_qc(mt).cols().flatten().key_by('s')
sample_qc_ht.export(output=SAMPLE_QC_OUT_FAILED)



