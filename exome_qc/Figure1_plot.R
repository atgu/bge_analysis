# ## 

library(ggplot2)
library(ggsci)
library(ggpubr)
library(cowplot)

style = read.csv('style_sheet_combinedGPC.csv', header=T)
dat_all = read.delim("pumas_bge_wave1_sample_qc_all_outputs_01-19-2024.tsv", header=T, sep="\t")
dat_all = merge(dat_all, style, by='TERRA_WORKSPACE')

dat_pass = subset(dat_all, sample_qc_outcome=="pass")

# exome cdf
exome_cov = ggplot(dat_pass, aes(exome_target_fraction_10x, y=1-..y.., color=Name, linetype=collection_method)) + 
	stat_ecdf(geom = "step", linewidth=1.2) + #scale_color_npg() + 
	theme_bw() + 
	theme(axis.text = element_text(size=14), axis.title = element_text(size=14), legend.text = element_text(size=14), legend.title=element_text(size=14)) +
	xlab("Fraction of Exome Target >=10x") +
	xlim(0.8, 1) + 
	ylab("Proportion of Samples") +
	# labs(color="Site") + 
	guides(linetype='none')  +
	guides(color='none')  +
	scale_linetype_manual(values=c("dotted","solid")) +
	 scale_color_manual(name = 'Cohort', 
		values = c("NeuroGAP_AddisEthiopia" = "#000075", 
					"NeuroGAP_CapeTownSouthAfrica" = "#3cb44b", 
					"NeuroGAP_MoiKenya" = "#42d4f4", 
					"NeuroGAP_KEMRIKenya" = "#4363d8", 
					"NeuroGAP_MakerereUganda" = "#aaffc3", 
					"GPC_USA" = "#e6194B", 
					"PaisaColombia" = "#f58231"))


## PCA plot
PCS = "/Users/juliasealock/Desktop/bge_manuscript/ancestry/agvp_awigen_ref/recombine_pca/recombine_datasets_pca__pc_scores.txt"
SAMPLES = '/Users/juliasealock/Desktop/bge_manuscript/ancestry/agvp_awigen_ref/recombine_pca/recombine_datasets_pca__PCA_samples.txt'
BGE_META = "bge_wave1_pumas_samples_manifest_11-17-2023.tsv"

pcs = read.delim(PCS)
samples = read.delim(SAMPLES, header=T)
bge_meta = read.delim(BGE_META, header=T)

pcs$scores = gsub("\\[|\\]", "", pcs$scores)
pcs = separate(data = pcs, col = scores, into = c("PC1","PC2","PC3","PC4","PC5","PC6","PC7","PC8","PC9","PC10"), sep = ",")
pcs[,2:11] <- sapply(pcs[,2:11],as.numeric)

bge_samples = subset(samples, Dataset=="PUMAS")
bge_pcs = pcs[(pcs$s %in% bge_samples$s),]
cohort = dat_pass[c(1,2,61:64)]
bge_pcs_cohort = merge(cohort, bge_pcs, by="s")


pc1_pc2 = ggplot() + 
     geom_point(data=bge_pcs_cohort, aes(x=PC1, y=PC2, color=Name, shape=Name))+
     theme_bw() + 
	theme(axis.text = element_text(size=14), axis.title = element_text(size=14), legend.text = element_text(size=14), legend.title = element_text(size=14)) +
    #  guides(color='none') + 	
	 scale_color_manual(name = 'Cohort', 
		values = c("NeuroGAP_AddisEthiopia" = "#000075", 
					"NeuroGAP_CapeTownSouthAfrica" = "#3cb44b", 
					"NeuroGAP_MoiKenya" = "#42d4f4", 
					"NeuroGAP_KEMRIKenya" = "#4363d8", 
					"NeuroGAP_MakerereUganda" = "#aaffc3", 
					"GPC_USA" = "#e6194B", 
					"PaisaColombia" = "#f58231")) +
	scale_shape_manual(name='Cohort',
		values = c("NeuroGAP_AddisEthiopia" = 16, 
					"NeuroGAP_CapeTownSouthAfrica" = 16, 
					"NeuroGAP_MoiKenya" = 16, 
					"NeuroGAP_KEMRIKenya" = 16, 
					"NeuroGAP_MakerereUganda" = 16, 
					"GPC_USA" = 17, 
					"PaisaColombia" = 15)) +
	guides(color = guide_legend(override.aes = list(size = 3))) 


pc2_pc3 = ggplot() + 
     geom_point(data=bge_pcs_cohort, aes(x=PC2, y=PC3, color=Name, shape=Name))+
     theme_bw() + 
     theme(axis.text = element_text(size=14), axis.title = element_text(size=14)) +
     guides(color='none', shape='none') + 	
	 scale_color_manual(name = 'Cohort', 
		values = c("NeuroGAP_AddisEthiopia" = "#000075", 
					"NeuroGAP_CapeTownSouthAfrica" = "#3cb44b", 
					"NeuroGAP_MoiKenya" = "#42d4f4", 
					"NeuroGAP_KEMRIKenya" = "#4363d8", 
					"NeuroGAP_MakerereUganda" = "#aaffc3", 
					"GPC_USA" = "#e6194B", 
					"PaisaColombia" = "#f58231")) +
	scale_shape_manual(name='Cohort',
		values = c("NeuroGAP_AddisEthiopia" = 16, 
					"NeuroGAP_CapeTownSouthAfrica" = 16, 
					"NeuroGAP_MoiKenya" = 16, 
					"NeuroGAP_KEMRIKenya" = 16, 
					"NeuroGAP_MakerereUganda" = 16, 
					"GPC_USA" = 17, 
					"PaisaColombia" = 15))

pc3_pc4 = ggplot() + 
     geom_point(data=bge_pcs_cohort, aes(x=PC3, y=PC4, color=Name, shape=Name))+
     theme_bw() + 
     theme(axis.text = element_text(size=14), axis.title = element_text(size=14)) +
    #  guides(color='none', shape='none') + 	
	 scale_color_manual(name = 'Cohort', 
		values = c("NeuroGAP_AddisEthiopia" = "#000075", 
					"NeuroGAP_CapeTownSouthAfrica" = "#3cb44b", 
					"NeuroGAP_MoiKenya" = "#42d4f4", 
					"NeuroGAP_KEMRIKenya" = "#4363d8", 
					"NeuroGAP_MakerereUganda" = "#aaffc3", 
					"GPC_USA" = "#e6194B", 
					"PaisaColombia" = "#f58231")) +
	scale_shape_manual(name='Cohort',
		values = c("NeuroGAP_AddisEthiopia" = 16, 
					"NeuroGAP_CapeTownSouthAfrica" = 16, 
					"NeuroGAP_MoiKenya" = 16, 
					"NeuroGAP_KEMRIKenya" = 16, 
					"NeuroGAP_MakerereUganda" = 16, 
					"GPC_USA" = 17, 
					"PaisaColombia" = 15))

	 




## sample qc metrics
### call rate, mean dp, mean gq, wgs depth 
dat_long = gather(dat_pass, metric, value, exome_target_fraction_10x:sample_qc.r_insertion_deletion)
##
dat2 = subset(dat_long, metric=="wgs_est_coverage" | metric=="sample_qc.call_rate" | metric=="sample_qc.dp_stats.mean" | metric=="sample_qc.gq_stats.mean")
dat2$metric = ifelse(dat2$metric=="sample_qc.call_rate", "Coding Call Rate", 
				ifelse(dat2$metric=="sample_qc.dp_stats.mean", "Mean Coding Depth", 
					ifelse(dat2$metric=="sample_qc.gq_stats.mean", "Mean Coding Genotype Quality", 
						ifelse(dat2$metric=="wgs_est_coverage", "Estimated WGS Coverage", "NA"))))

dat2$metric = factor(dat2$metric, levels=c("Estimated WGS Coverage", "Mean Coding Depth", "Coding Call Rate", "Mean Coding Genotype Quality"))
sample_qc = ggplot(dat2, aes(x=Name, y=value)) +
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
	    facet_wrap(~metric, ncol=4, scales="free_x") +
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
					"PaisaColombia" = "#f58231")) 



## combine into one plot

my_plot = ggarrange(pc1_pc2, pc3_pc4, exome_cov, ncol=3, common.legend = TRUE, labels = c("A","B","C"))
out = plot_grid(my_plot, sample_qc, labels = c('', 'D'), nrow=2)


png(paste0("pumas_figure1.png"), width = 18, height = 12, units = 'in', res = 300)
print(out)
dev.off()


pdf(paste0("pumas_figure1.pdf"), width = 18, height = 12)
print(out)
dev.off()

