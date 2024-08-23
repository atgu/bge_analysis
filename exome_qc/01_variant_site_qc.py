## exome QC
import hail as hl
hl.init(driver_cores=8, worker_memory='highmem', tmp_dir="gs://schema_jsealock/tmp/")

VCF_PATH = 'BGE_Exome_Callset_Global_NovaSeq6000s_chr*.*.hard_filtered_with_genotypes.vcf.gz'
TARGET_INTERVALS = "Twist_Alliance_Clinical_Research_Exome_Covered_Targets_hg38-34.9MB.bed" ## no padding
LCR_PATH = "LCRFromHengHg38.bed"
VQSR_VCF = 'BGE_Exome_Callset_Global_NovaSeq6000s.filtered.*.vcf.gz'
MT_OUT = 'bge_wave1_variant_qc.mt'

## 
print("read in vcfs")
vcfs = [entry['path'] for entry in hl.hadoop_ls(VCF_PATH)]
mt = hl.import_vcf(vcfs, reference_genome='GRCh38', force_bgz=True, array_elements_required=False)

print("filter multi alleles and split")
mt = mt.filter_rows(mt.alleles.length() <= 6)

bi = mt.filter_rows(hl.len(mt.alleles) == 2)
bi = bi.annotate_rows(a_index=1, was_split=False)
multi = mt.filter_rows(hl.len(mt.alleles) > 2)
split = hl.split_multi_hts(multi)
mt = split.union_rows(bi)

print("pre filter variants")
mt = mt.filter_entries(
    hl.is_defined(mt.GT) &
    (
        (mt.GT.is_hom_ref() & 
            (
                # ((mt.AD[0] / mt.DP) < 0.8) | # Has to be removed because allele depth no longer defined for hom ref calls.
                (mt.GQ < 20) |
                (mt.DP < 10)
        	)
        ) |
        (mt.GT.is_het() & 
        	( 
                # (((mt.AD[0] + mt.AD[1]) / mt.DP) < 0.8) |  
                ((mt.AD[1] / mt.DP) < 0.2) | 
                ((mt.AD[1] / mt.DP) > 0.8) | 
                (mt.PL[0] < 20) |
                (mt.DP < 10)
        	)
        ) |
        (mt.GT.is_hom_var() & 
        	(
                ((mt.AD[1] / mt.DP) < 0.8) |
                (mt.PL[0] < 20) |
                (mt.DP < 10)
        	)
        )
    ),
    keep = False
)


intervals = [hl.parse_locus_interval(x, reference_genome='GRCh38') for x in ['chr1:START-chr22:END', 'chrX:START-chrX:END', 'chrY:START-chrY:END']]
mt = hl.filter_intervals(mt, intervals)


print("vqsr filter")
part2mt = hl.import_vcf(VQSR_VCF, reference_genome='GRCh38', force_bgz=True)
part2mt = part2mt.annotate_rows(fail_VQSR = hl.len(part2mt.filters) != 0)
passed = part2mt.filter_rows(part2mt.fail_VQSR, keep=False).rows()
mt = mt.filter_rows(hl.is_defined(passed[mt.row_key]), keep=True)

#filter out LCRs 
print("lcr filter")
lcr_intervals = hl.import_locus_intervals(LCR_PATH, reference_genome='GRCh38', skip_invalid_intervals=True)
mt = mt.filter_rows(hl.is_defined(lcr_intervals[mt.locus]), keep=False)

# interval filter
print("interval filter")
intervals = hl.import_locus_intervals(TARGET_INTERVALS, reference_genome="GRCh38")
mt = mt.filter_rows(hl.is_defined(intervals[mt.locus]), keep=True)

# Filter out the invariant rows.
mt = hl.variant_qc(mt, name='qc')
mt = mt.filter_rows((mt.qc.AF[0] > 0.0) & (mt.qc.AF[0] < 1.0))

# make into minimal representation 
print("min rep")
mt = mt.annotate_rows(min_rep = hl.min_rep(mt.locus, mt.alleles))
mt = mt.key_rows_by('min_rep')
mt = mt.drop('locus', 'alleles')

mt = mt.annotate_rows(locus = mt.min_rep.locus,
                      alleles = mt.min_rep.alleles)
mt = mt.key_rows_by('locus', 'alleles')
mt = mt.drop('min_rep')


print("write prefiltered mt")
mt.write(MT_OUT, overwrite=True)


