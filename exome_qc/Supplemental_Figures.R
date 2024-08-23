## supplemental plots for pcs + ref panel and sample qc metrics

library(ggplot2)
library(ggsci)
library(ggpubr)
library(cowplot)

style = read.csv('/Users/juliasealock/Desktop/bge_manuscript/style_sheet_combinedGPC.csv', header=T)
dat_all = read.delim("/Users/juliasealock/Desktop/bge_manuscript/pumas_bge_wave1_sample_qc_all_outputs_01-19-2024.tsv", header=T, sep="\t")
dat_all = merge(dat_all, style, by='TERRA_WORKSPACE')

dat_pass = subset(dat_all, sample_qc_outcome=="pass")

## PCA plot

PCS = "/Users/juliasealock/Desktop/bge_manuscript/ancestry/agvp_awigen_ref/recombine_pca/recombine_datasets_pca__pc_scores.txt"
SAMPLES = '/Users/juliasealock/Desktop/bge_manuscript/ancestry/agvp_awigen_ref/recombine_pca/recombine_datasets_pca__PCA_samples.txt'
BGE_META = "/Users/juliasealock/Desktop/bge_manuscript/bge_wave1_pumas_samples_manifest_11-17-2023.tsv"
META = '/Users/juliasealock/Desktop/bge_manuscript/ancestry/agvp_awigen_ref/global_pca_metadata.tsv'


samples = read.delim(SAMPLES, header=T)
bge_meta = read.delim(BGE_META, header=T)


bge_samples = subset(samples, Dataset=="PUMAS")
bge_pcs = pcs[(pcs$s %in% bge_samples$s),]
cohort = dat_pass[c(1,2,61:64)]
bge_pcs_cohort = merge(cohort, bge_pcs, by="s")


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
ref_pcs_meta = merge(ref_pcs, meta, by.x='s', by.y='IID')
ref_pcs_meta$continent = ifelse(ref_pcs_meta$Region=="East Africa","AFR",
                          ifelse(ref_pcs_meta$Region=="Southern Africa", "AFR", 
                            ifelse(ref_pcs_meta$Region=="West Africa", "AFR", ref_pcs_meta$Region)))

bge_pcs_cohort$Name = ifelse(bge_pcs_cohort$Name=='PaisaColombia', 'Paisa_Colombia', bge_pcs_cohort$Name)
bge_pcs_cohort$Name = factor(bge_pcs_cohort$Name, levels=c('NeuroGAP_AddisEthiopia','NeuroGAP_KEMRIKenya','NeuroGAP_MoiKenya','NeuroGAP_CapeTownSouthAfrica','NeuroGAP_MakerereUganda','Paisa_Colombia','GPC_USA'))


# png(paste0("pumas_supplemental_pc1_pc2_with_ref_all_samples.png"), width = 12, height = 12, units = 'in', res = 300)
# pc1_pc2 = ggplot() + 
#     geom_point(data=ref_pcs_meta, aes(x=PC1, y=PC2, fill=continent), alpha=0.7, shape=24, stroke=NA, size=3) +
#     geom_point(data=bge_pcs_cohort, aes(x=PC1, y=PC2, color=Name), alpha=0.7) +
#     theme_bw() + 
# 	theme(axis.text = element_text(size=14), axis.title = element_text(size=14), legend.text = element_text(size=14), legend.title = element_text(size=14)) +
# 	 scale_color_manual(name = 'Cohort', 
# 		values = c("NeuroGAP_AddisEthiopia" = "#000075", 
# 					"NeuroGAP_CapeTownSouthAfrica" = "#3cb44b", 
# 					"NeuroGAP_MoiKenya" = "#42d4f4", 
# 					"NeuroGAP_KEMRIKenya" = "#4363d8", 
# 					"NeuroGAP_MakerereUganda" = "#aaffc3", 
# 					"GPC_USA" = "#e6194B", 
# 					"PaisaColombia" = "#f58231")) +
#     scale_fill_manual(name='Reference Samples',
#         values = c('AFR' = '#ff1493', 
#                     'AMR' = '#00ffff',
#                     'CSA' = '#ffff00',
#                     'EAS' = '#8b0000',
#                     'EUR' = '#00ff00',
#                     'MID' = '#ff4500', 
#                     'OCE' = '#800080')) +
# 	guides(color = guide_legend(override.aes = list(size = 3))) 
# print(pc1_pc2)
# dev.off()

# bge_pcs_cohort$Name = factor(bge_pcs_cohort$Name, levels=c("NeuroGAP_AddisEthiopia","NeuroGAP_CapeTownSouthAfrica","NeuroGAP_MoiKenya","NeuroGAP_KEMRIKenya","NeuroGAP_MakerereUganda","GPC_USA","PaisaColombia"))

# png(paste0("pumas_supplemental_pc1_pc2_with_ref_continent_by_cohort.png"), width = 15, height = 12, units = 'in', res = 300)
# pc1_pc2 = ggplot() + 
#     geom_point(data=ref_pcs_meta, aes(x=PC1, y=PC2, fill=continent), shape=24, stroke=NA, size=2, alpha=0.7) +
#     geom_point(data=bge_pcs_cohort, aes(x=PC1, y=PC2), color='gold', alpha=0.7) +
#     theme_bw() + 
# 	theme(axis.text = element_text(size=14), axis.title = element_text(size=14), legend.text = element_text(size=14), legend.title = element_text(size=14)) +
#     scale_fill_d3() + 
#     facet_wrap(~Name) + 
#     guides(fill=guide_legend('Reference Sample Continent')) + 
#     guides(fill = guide_legend(override.aes = list(size = 4))) 
# print(pc1_pc2)
# dev.off()


# png(paste0("pumas_supplemental_pc1_pc2_with_ref_region_by_cohort.png"), width = 15, height = 12, units = 'in', res = 300)
pc1_pc2 = ggplot() + 
    geom_point(data=ref_pcs_meta, aes(x=PC1, y=PC2, fill=Region), shape=24, stroke=NA, size=2, alpha=0.7) +
    geom_point(data=bge_pcs_cohort, aes(x=PC1, y=PC2), color='gold', alpha=0.7) +
    theme_bw() + 
	theme(axis.text = element_text(size=14), axis.title = element_text(size=14), legend.text = element_text(size=14), legend.title = element_text(size=14)) +
    scale_fill_d3() + 
    facet_wrap(~Name, ncol=4) + 
    guides(fill=guide_legend('Reference Sample Region')) + 
    guides(fill = guide_legend(override.aes = list(size = 4))) +
		    theme(
	    	panel.grid.major = element_blank(),
        # 	# panel.grid.minor = element_blank(),
        	strip.background=element_rect(fill="white"),
        	panel.border = element_rect(colour = "black", fill = NA)) +
		theme(legend.position = "top") +
		theme(strip.text = element_text(size = 14))
# print(pc1_pc2)
# dev.off()



# NeuroGAP with all of the reference panels, and then Paisa and GPC just with HGDP1KG, to make it easier to see the non-AFR dots
ng = subset(bge_pcs_cohort, Name!='GPC_USA' & Name!='Paisa_Colombia')
gpc_paisa = subset(bge_pcs_cohort, Name=='GPC_USA' | Name=='Paisa_Colombia')


pc1_pc2_ng = ggplot() + 
    geom_point(data=ref_pcs_meta, aes(x=PC1, y=PC2, fill=Region), shape=24, stroke=NA, size=2, alpha=0.7) +
    geom_point(data=ng, aes(x=PC1, y=PC2), color='gold', alpha=0.7) +
    theme_bw() + 
	theme(axis.text = element_text(size=14), axis.title = element_text(size=14), legend.text = element_text(size=14), legend.title = element_text(size=14)) +
    scale_fill_d3() + 
    facet_wrap(~Name, ncol=3) + 
    guides(fill=guide_legend('Reference Sample Region')) + 
    guides(fill = guide_legend(override.aes = list(size = 4))) +
		    theme(
	    	panel.grid.major = element_blank(),
        # 	# panel.grid.minor = element_blank(),
        	strip.background=element_rect(fill="white"),
        	panel.border = element_rect(colour = "black", fill = NA)) +
		theme(legend.position = "top") +
		theme(strip.text = element_text(size = 14))

pc3_pc4_ng = ggplot() + 
    geom_point(data=ref_pcs_meta, aes(x=PC3, y=PC4, fill=Region), shape=24, stroke=NA, size=2, alpha=0.7) +
    geom_point(data=ng, aes(x=PC1, y=PC2), color='gold', alpha=0.7) +
    theme_bw() + 
	theme(axis.text = element_text(size=14), axis.title = element_text(size=14), legend.text = element_text(size=14), legend.title = element_text(size=14)) +
    scale_fill_d3() + 
    facet_wrap(~Name, ncol=3) + 
    guides(fill=guide_legend('Reference Sample Region')) + 
    guides(fill = guide_legend(override.aes = list(size = 4))) +
		    theme(
	    	panel.grid.major = element_blank(),
        # 	# panel.grid.minor = element_blank(),
        	strip.background=element_rect(fill="white"),
        	panel.border = element_rect(colour = "black", fill = NA)) +
		theme(legend.position = "top") +
		theme(strip.text = element_text(size = 14))

ref_pcs_meta2 = subset(ref_pcs_meta, source=='1000 Genomes Project' | source=='HGDP')
pc1_pc2_paisa = ggplot() + 
    geom_point(data=ref_pcs_meta2, aes(x=PC1, y=PC2, fill=Region), shape=24, stroke=NA, size=2, alpha=0.7) +
    geom_point(data=gpc_paisa, aes(x=PC1, y=PC2), color='gold', alpha=0.7) +
    theme_bw() + 
	theme(axis.text = element_text(size=14), axis.title = element_text(size=14), legend.text = element_text(size=14), legend.title = element_text(size=14)) +
    scale_fill_d3() + 
    facet_wrap(~Name, ncol=1) + 
    guides(fill=guide_legend('Reference Sample Region')) + 
    guides(fill = guide_legend(override.aes = list(size = 4))) +
		    theme(
	    	panel.grid.major = element_blank(),
        # 	# panel.grid.minor = element_blank(),
        	strip.background=element_rect(fill="white"),
        	panel.border = element_rect(colour = "black", fill = NA)) +
		theme(legend.position = "top") +
		theme(strip.text = element_text(size = 14))

pc3_pc4_paisa = ggplot() + 
    geom_point(data=ref_pcs_meta2, aes(x=PC3, y=PC4, fill=Region), shape=24, stroke=NA, size=2, alpha=0.7) +
    geom_point(data=gpc_paisa, aes(x=PC1, y=PC2), color='gold', alpha=0.7) +
    theme_bw() + 
	theme(axis.text = element_text(size=14), axis.title = element_text(size=14), legend.text = element_text(size=14), legend.title = element_text(size=14)) +
    scale_fill_d3() + 
    facet_wrap(~Name, ncol=1) + 
    guides(fill=guide_legend('Reference Sample Region')) + 
    guides(fill = guide_legend(override.aes = list(size = 4))) +
		    theme(
	    	panel.grid.major = element_blank(),
        # 	# panel.grid.minor = element_blank(),
        	strip.background=element_rect(fill="white"),
        	panel.border = element_rect(colour = "black", fill = NA)) +
		theme(legend.position = "top") +
		theme(strip.text = element_text(size = 14))




## sample qc metrics 
dat_long = gather(dat_pass, metric, value, chimeras:sample_qc.r_insertion_deletion)
##
metrics = c("chimeras", 'contamination','sample_qc.n_deletion','sample_qc.n_insertion','sample_qc.n_singleton', 'sample_qc.n_transition','sample_qc.n_transversion','sample_qc.r_het_hom_var','sample_qc.r_insertion_deletion','sample_qc.r_ti_tv')
dat2 = dat_long[(dat_long$metric %in% metrics),]

names = c('Chimeric Read Rate', 'Contamination', 'N Deletions', 'N Insertions', 'N Singletons', 'N Transitions', 'N Transversions','Heterozygosity Ratio','Indel Ratio', 'Ti/Tv Ratio')
names = data.frame(metrics, names)

dat3 = merge(dat2, names, by.x='metric', by.y='metrics')
dat3$Name = ifelse(dat3$Name=='PaisaColombia', 'Paisa_Colombia', dat3$Name)
dat3$Name = factor(dat3$Name, levels=c('NeuroGAP_AddisEthiopia','NeuroGAP_KEMRIKenya','NeuroGAP_MoiKenya','NeuroGAP_CapeTownSouthAfrica','NeuroGAP_MakerereUganda','Paisa_Colombia','GPC_USA'))


# png(paste0("pumas_supplemental_sample_qc_passing_samples.png"), width = 15, height = 10, units = 'in', res = 300)
sample_qc = ggplot(dat3, aes(x=Name, y=value)) +
	    geom_boxplot(outlier.shape = NA) + 
	    geom_jitter(width=0.2, height=0, alpha=0.5, aes(color=Name)) + 
	    coord_flip() +
	    theme_bw() +
	    theme(axis.title.x = element_text(margin = ggplot2::margin(t=10))) +
	    # labs(color="Ancestry") +
		guides(color='none') + 
	    ylab("") +
	    xlab(" ") +
	    theme(axis.text=element_text(size=12),
	        axis.title=element_text(size=12)) +
	    scale_color_d3() +
	    facet_wrap(~names, ncol=5, scales="free_x") +
	    scale_x_discrete(limits=rev) +
	    theme(strip.text = element_text(size = 14)) +
	    theme(
	    	panel.grid.major = element_blank(),
        # 	# panel.grid.minor = element_blank(),
        	strip.background=element_rect(fill="white"),
        	panel.border = element_rect(colour = "black", fill = NA)) +
		scale_color_manual(name = 'Cohort', 
			values = c("NeuroGAP_AddisEthiopia" = "#000075", 
					"NeuroGAP_CapeTownSouthAfrica" = "#3cb44b", 
					"NeuroGAP_MoiKenya" = "#42d4f4", 
					"NeuroGAP_KEMRIKenya" = "#4363d8", 
					"NeuroGAP_MakerereUganda" = "#aaffc3", 
					"GPC_USA" = "#e6194B", 
					"Paisa_Colombia" = "#f58231")) 
print(sample_qc)
dev.off()


pc1_pc2_hex = ggplot(data=bge_pcs_cohort, aes(x=PC1, y=PC2)) + 
    theme_bw() + 
	theme(axis.text = element_text(size=14), axis.title = element_text(size=14), legend.text = element_text(size=14), legend.title = element_text(size=14)) +
    geom_hex(bins=75) +
    scale_fill_gradientn(name='Count',colours = rev(brewer.pal(5,'Spectral'))) +
	theme(plot.margin = margin(0, 2, 0, 2, "cm"))
	# facet_wrap(~Name)

pc3_pc4_hex = ggplot(data=bge_pcs_cohort, aes(x=PC3, y=PC4)) + 
    theme_bw() + 
	theme(axis.text = element_text(size=14), axis.title = element_text(size=14), legend.text = element_text(size=14), legend.title = element_text(size=14)) +
    geom_hex(bins=75) +
    scale_fill_gradientn(name='Count',colours = rev(brewer.pal(5,'Spectral'))) +
	theme(plot.margin = margin(0, 2, 0, 2, "cm"))
	# facet_wrap(~Name)




# pc1_plots = ggarrange(pc1_pc2_ng, pc1_pc2_paisa, ncol=2, common.legend=T, labels=c('A.','C.'), widths=c(1.75, 0.8))
# pc3_plots = ggarrange(pc3_pc4_ng, pc3_pc4_paisa, ncol=2, common.legend=T, labels=c('B.','D.'), widths=c(1.75, 0.8))

# pc_plots = ggarrange(pc1_plots, pc3_plots, common.legend=TRUE, ncol=1)

pc_plots = ggarrange(pc1_pc2_ng, pc1_pc2_paisa, pc3_pc4_ng, pc3_pc4_paisa, common.legend=TRUE, ncol=2, nrow=2, widths=c(2, 0.75), labels=c('A.','','B.',''))

p2 = ggarrange(pc_plots, sample_qc, ncol=1, heights=c(1, 0.75), labels=c('', 'C.'))

png(paste0("pumas_supplemental_sample_qc_figure.png"), width = 15, height = 18, units = 'in', res = 300)
print(p2)
dev.off()


hex = ggarrange(pc1_pc2_hex, pc3_pc4_hex, common.legend=T, legend='right', widths = c(0.5, 0.5), heights = c(0.5, 0.5), labels=c('C.','D.'))



p2 = ggarrange(pc_plots, hex, sample_qc, ncol=1, heights=c(1, 0.33, 0.75), labels=c('', '', 'E.'))

png(paste0("pumas_supplemental_sample_qc_figure.png"), width = 15, height = 20, units = 'in', res = 300)
print(p2)
dev.off()

