## merge AWIGEN + AGVP + 1KGP + HGDP + PUMAS
## filter and ld prune
## run pca

import hail as hl 
PCA_OUT = 'recombine_datasets_pca_'

# lift over agvp
agvp = hl.import_plink(bed='AGVP.post.autosomes.unrelated.bed', bim='AGVP.post.autosomes.unrelated.bim', fam='AGVP.post.autosomes.unrelated.fam', reference_genome="GRCh37")
rg37 = hl.get_reference('GRCh37')  
rg38 = hl.get_reference('GRCh38')  
rg37.add_liftover('gs://hail-common/references/grch37_to_grch38.over.chain.gz', rg38)  
agvp = agvp.annotate_rows(new_locus=hl.liftover(agvp.locus, 'GRCh38', include_strand=True))#, old_locus=ht.locus)  
agvp = agvp.filter_rows(hl.is_defined(agvp.new_locus) & ~agvp.new_locus.is_negative_strand)  
agvp = agvp.key_rows_by(locus=agvp.new_locus.result, alleles=agvp.alleles)  
agvp = agvp.annotate_cols(Dataset="AGVP")

# load awigen
awigen = hl.import_plink(bed="awigen-annot-lifted.bed", bim="awigen-annot-lifted.bim",fam="awigen-annot-lifted.fam",reference_genome='GRCh38')
awigen = awigen.annotate_cols(Dataset="AWIGEN")

# load 1kgp + hdgp
gnomad1 = hl.experimental.load_dataset(name='gnomad_hgdp_1kg_subset_dense',
                                  version='3.1.2',
                                  reference_genome='GRCh38',
                                  region='us',
                                  cloud='gcp') # https://hail.is/docs/0.2/experimental/index.html#hail.experimental.load_dataset
gnomad = gnomad1.annotate_cols(Dataset = 'gnomad_hgdp_1kg')
# Take subset of entry fields 
gnomad = gnomad.select_entries('DP', 'GQ', 'GT')
# Take only sites passing VQSR
gnomad = gnomad.filter_rows(gnomad.filters == hl.empty_set(hl.tstr), keep = True)

# load pumas
bge = hl.read_matrix_table('bge_wave1_variant_qc.mt')
meta_bge = hl.import_table('bge_wave1_manuscript_samples.tsv', key="SAMPLE_ALIAS")
bge = bge.filter_cols(hl.is_defined(meta_bge[bge.col_key]))
samples = hl.import_table('bge_manuscript_filtered_samples_chim_contam_5_exome_90_wgs_1x.tsv', key='s')
bge = bge.filter_cols(hl.is_defined(samples[bge.col_key]))
bge = bge.annotate_cols(Dataset="PUMAS")

## select cols, rows, GT field
agvp = agvp.select_cols('Dataset').select_rows().select_entries("GT")
awigen = awigen.select_cols('Dataset').select_rows().select_entries("GT")
gnomad = gnomad.select_cols('Dataset').select_rows().select_entries("GT").select_globals()
bge = bge.select_cols('Dataset').select_rows().select_entries("GT").select_globals()

# merge datasets
agvp_awigen = agvp.union_cols(awigen)
agvp_awigen_gnomad = agvp_awigen.union_cols(gnomad)
mt = agvp_awigen_gnomad.union_cols(bge)

## call rate > 0.99, MAF > 0.1%
mt = hl.variant_qc(mt)
mt = mt.filter_rows(mt.variant_qc.call_rate>0.99)

mt = mt.filter_rows(hl.len(mt.alleles) == 2)
mt = mt.filter_rows(((mt.variant_qc.AF[0] > 0.001) & (mt.variant_qc.AF[1] > 0.001)) & ((mt.variant_qc.AF[0] < 0.999) & (mt.variant_qc.AF[1] < 0.999)))

print("checkpoint 1")
mt = mt.checkpoint(f'{PCA_OUT}.mt', overwrite=True)

mt = hl.read_matrix_table(f'{PCA_OUT}.mt')
print(mt.count())

# Remove ChrX to avoid errors 
mt = mt.filter_rows(mt.locus.contig != "chrX")

print("ld prune")
pruned_variant_table = hl.ld_prune(mt.GT, r2=0.1, bp_window_size=500000)
pruned_variant_table.write(f'{PCA_OUT}pruned_variant_table.ht')

mt = mt.filter_rows(hl.is_defined(pruned_variant_table[mt.row_key]))

mt = mt.checkpoint(f'{PCA_OUT}_ldpruned.mt', overwrite=True)

print("save samples lists")
mt.cols().export(f'{PCA_OUT}_PCA_samples.txt') 

print("run pca")
# # run pca
_, score_table, _ = hl.hwe_normalized_pca(mt.GT, k=10, compute_loadings=False)
score_table.flatten().export(f'{PCA_OUT}_pc_scores.txt')
