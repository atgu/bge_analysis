
library(tidyr)
library(tidyverse)
library(randomForest)
library(ggsci)
library(RColorBrewer)
library(colorRamps)

## recombined pca scores

PCS = 'recombine_datasets_pca__pc_scores.txt'
SAMPLES = 'recombine_datasets_pca__PCA_samples.txt'
META = 'global_pca_metadata.tsv'
BGE_META = "bge_wave1_pumas_samples_manifest_11-17-2023.tsv"

pcs = read.delim(PCS)
samples = read.delim(SAMPLES, header=T)
meta = read.delim(META, header=T)
meta = subset(meta, source!="NeuroGAP") # remove neurogap from ref 
bge_meta = read.delim(BGE_META, header=T)

pcs$scores = gsub("\\[|\\]", "", pcs$scores)
pcs = separate(data = pcs, col = scores, into = c("PC1","PC2","PC3","PC4","PC5","PC6","PC7","PC8","PC9","PC10"), sep = ",")
pcs[,2:11] <- sapply(pcs[,2:11],as.numeric)

ref_samples = subset(samples, Dataset=="AGVP" | Dataset=="AWIGEN" | Dataset=="gnomad_hgdp_1kg")
bge_samples = subset(samples, Dataset=="PUMAS")

ref_pcs = pcs[(pcs$s %in% ref_samples$s),]
bge_pcs = pcs[(pcs$s %in% bge_samples$s),]

### collapse all AFR regions into 1
##### These are the results I used for ancestry determination 
ref_pcs_meta$continent = ifelse(ref_pcs_meta$Region=="East Africa","AFR",
                          ifelse(ref_pcs_meta$Region=="Southern Africa", "AFR", 
                            ifelse(ref_pcs_meta$Region=="West Africa", "AFR", ref_pcs_meta$Region)))



forest = randomForest(as.factor(continent) ~ PC1 + PC2 + PC3 + PC4 + PC5 + PC6 + PC7 + PC8 + PC9 + PC10,
                         data = subset(ref_pcs_meta, 1 == 1),
                         importance = TRUE,
                         # ntree = 10000)
                         ntree = 1000)

# Do the prediction
fit_data_all = data.frame(predict(forest, bge_pcs, type='prob'), Sample = bge_pcs$s)

# Do some re-organization voodoo
fit_data = fit_data_all %>%
  gather(Population, Probability, -Sample) %>%
  group_by(Sample) %>%
  slice(which.max(Probability))

fit_cohort = merge(fit_data, tw, by.x="Sample", by.y="SAMPLE_ALIAS")
bge_pcs_fit = merge(bge_pcs, fit_cohort, by.x='s', by.y="Sample")


write.table(bge_pcs_fit, "pumas_bge_samples_pcs_rf_prob.txt", col.names=T, row.names=F, quote=F, sep="\t")
