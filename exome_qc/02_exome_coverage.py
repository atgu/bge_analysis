

import hail as hl

hl.init(driver_cores=8, worker_memory='highmem', tmp_dir="gs://schema_jsealock/tmp/")

VCF_PATH = 'BGE_Exome_Callset_Global_NovaSeq6000s_chr*.*.hard_filtered_with_genotypes.vcf.gz'
TARGET_INTERVALS = "Twist_Alliance_Clinical_Research_Exome_Covered_Targets_hg38-34.9MB.bed" ## no padding
META = 'bge_wave1_manuscript_samples.tsv'

## 
print("read in vcfs")
vcfs = [entry['path'] for entry in hl.hadoop_ls(VCF_PATH)]
mt = hl.import_vcf(vcfs, reference_genome='GRCh38', force_bgz=True, array_elements_required=False)

print("filter to pumas samples")
meta = hl.import_table(META, key="SAMPLE_ALIAS")
mt = mt.filter_cols(hl.is_defined(meta[mt.col_key]))

# interval filter
print("interval filter")
intervals = hl.import_locus_intervals(TARGET_INTERVALS, reference_genome="GRCh38")
mt = mt.filter_rows(hl.is_defined(intervals[mt.locus]), keep=True)

print("calc coverage")
ht = mt.annotate_cols(fraction_target_10x = hl.agg.fraction(mt.DP>=10)).cols()
ht.export('pumas_bge_wave1_fraction_twist_target_10x_no_qc.tsv.bgz')

