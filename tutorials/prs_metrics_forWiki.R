library(data.table)
library(pROC) #for calculating AUC & its 95% CIs
library(boot) #for calculating CIs for NKr2 & h2l_NKr2
library(ggplot2)

cal_NKr2_auc <- function(dat){
  ##base model with covariates only: Age + Sex + 10PCs
  glm0 <- glm(as.formula(paste("PHENO1 ~ ", exp0)), data = dat, family = binomial(logit)) 
  
  # logistic full model with PRS
  glm1 <- glm(as.formula(paste("PHENO1 ~ ", exp1)), data = dat, family = binomial(logit)) 
  
  # Calculate Cox & Snell R2 using log Likelihoods 
  LL1 <-  logLik(glm1)
  LL0 <-  logLik(glm0)
  CSr2 <-  round(1 - exp((2 / N) * (LL0[1] - LL1[1])), 6)
  
  # Calculate Nagelkerke's R2
  NKr2 <- round(CSr2 / (1 - exp((2 / N) * LL0[1])), 6)  
  # Test whether NKr2 is significantly different from 0
  devdiff <- round(glm0$deviance - glm1$deviance, 1) #Difference in deviance attributable to PRS
  df <- glm0$df.residual - glm1$df.residual #1 degree of freedom for single variable PRS
  NKr2_pval <- pchisq(devdiff, df, lower.tail = F)
  
  ## Calculate AUC using full model (PRS+covariates)
  auc1 <- round(auc(dat$PHENO1, glm1$linear.predictors), 3)
  
  glm2 <- glm(PHENO1 ~ ZSCORE, data = dat, family = binomial(logit))
  auc2 <- round(auc(dat$PHENO1, glm2$linear.predictors), 3)
  
  auc1_2.5 <- round(ci.auc(dat$PHENO1, glm1$linear.predictors)[1], 3)
  auc1_97.5 <- round(ci.auc(dat$PHENO1, glm1$linear.predictors)[3], 3)
  
  auc2_2.5 <- round(ci.auc(dat$PHENO1, glm2$linear.predictors)[1], 3)
  auc2_97.5 <- round(ci.auc(dat$PHENO1, glm2$linear.predictors)[3], 3)
  
  glm1.zbeta <- glm1$coefficients[2]
  glm2.zbeta <- glm2$coefficients[2]
  
  return(data.frame(NKr2, NKr2_pval,
                    auc1, auc1_2.5, auc1_97.5,
                    auc2, auc2_2.5, auc2_97.5, glm1.zbeta, glm2.zbeta))
}  

#################calculate proportion of variance explained on the liability scale#################
# Ref: Lee et al., Genet Epidemiol. 2012 Apr;36(3):214-24.
h2l_R2s <- function(k, r2, p) {
  # k baseline disease risk
  # r2 from a linear regression model of genomic profile risk score (R2v/NKr2)
  # p proportion of cases
  x <- qnorm(1 - k)
  z <- dnorm(x)
  i <- z / k
  C <- k * (1 - k) * k * (1 - k) / (z^2 * p * (1 - p))
  theta <- i * ((p - k) / (1 - k)) * (i * ((p - k) / (1 - k)) - x)
  e <- 1 - p^(2 * p) * (1 - p)^(2 * (1 - p))
  h2l_NKr2 <- C * e * r2 / (1 + C * e * theta * r2)
  h2l_NKr2 <- round(h2l_NKr2, 6)
  return(h2l_NKr2)
} 

# Calculate 95% CI for Nagelkerke's R2 & h2l_NKr2
boot_function <- function(data, indices){
  d <- data[indices,]
  NKr2 <- cal_NKr2_auc(d)$NKr2
  h2l_NKr2 <- h2l_R2s(K, NKr2, P)
  return(c(NKr2, h2l_NKr2))
}

pcvecs <- paste("PC",seq(1,10),sep="")
expnull <- "1"
exp0 <- paste0("SEX + AGE + ",paste0(pcvecs, collapse = " + "))
exp1 <- paste0("ZSCORE + SEX + AGE + ", (paste0(pcvecs, collapse = " + ")))

scorefile <- fread("AAU_scores.sscore")
head(scorefile)

pc_file <- fread("~/Desktop/bge_analysis/neuroGAP/neurogap_all_sites_autosomes_pcscores.txt")
pc_file$IID <- gsub(" ", "-",pc_file$s)
pc_file$IID <- gsub("_", "-",pc_file$IID)
pc_file$IID <- gsub("\\.", "-",pc_file$IID)
head(pc_file)

fam_file <- fread("~/Desktop/bge_analysis/PGS/neuroGAP_PGS/forWiki/plink/AAU_prs_qc2.fam")
colnames(fam_file) <- c("FID","IID","MID","PID","SEX","PHENO")
fam_file$MID <- NULL
fam_file$PID <- NULL
fam_file$SEX <- as.numeric(fam_file$SEX)
fam_file$PHENO <- as.numeric(fam_file$PHENO)
head(fam_file)

manifest.df <- fread("~/Desktop/bge_analysis/BGE_HailCallset_Wave2_Manifest_02202024.tsv")

#view the different project names
levels(as.factor(manifest.df$COHORT))

#subset to the project of interest (in our case "NeuroGAP-Psychosis_AAU")
manifest.sub <- subset(manifest.df, manifest.df$COHORT == "NeuroGAP-Psychosis_AAU")
manifest.sub$IID <- gsub(" ", "-",manifest.sub$SUBJECT_ID)
manifest.sub$IID <- gsub("_", "-",manifest.sub$IID)
manifest.sub$IID <- gsub("\\.", "-",manifest.sub$IID)

age.df <- as.data.frame(cbind(manifest.sub$IID, manifest.sub$AGE_UNSPECIFIED))
colnames(age.df) <- c("IID","AGE")
age.df$AGE <- as.numeric(age.df$AGE)
head(age.df)

#make PRS dataframe, remove unnecessary columns
prs <- merge(scorefile, fam_file, by="IID")
prs <- merge(prs,pc_file,by="IID")
prs <- merge(prs,age.df,by="IID")
prs$`#FID` <- NULL

#PHENO uses plink coding 1=control, 2=case
#lets change PHENO1 to be 0=control, 1=case
prs$PHENO1 <- prs$PHENO - 1

#get zscore (centered by controls)
scores <- c("SCORESUM", "SCORE1_SUM")
prs[,SCORESUM := get(grep(paste(scores, collapse = "|"), names(prs), value = T))]
prs[,ZSCORE := (SCORESUM - mean(SCORESUM[PHENO1 == 0]))/sd(SCORESUM[PHENO1==0])] # z-score of the profile score

N <-  prs[PHENO1 %in% c(0,1), .N]

tmp <- cal_NKr2_auc(prs)
NKr2 <- tmp$NKr2
NKr2_pval <- tmp$NKr2_pval

K=0.01 #prevalence of SCZ 

P <- prs[PHENO1 == 1, .N] / N ##proportion of cases

h2l_NKr2 <- h2l_R2s(K, NKr2, P)

results <- boot(prs, boot_function, R = 1000)

NKr2_2.5 <- round(boot.ci(results, type ="perc", index = 1)[[4]][1,][4], 6)
NKr2_97.5 <- round(boot.ci(results, type ="perc", index = 1)[[4]][1,][5], 6)
h2l_NKr2_2.5 <- round(boot.ci(results, type ="perc", index = 2)[[4]][1,][4], 6)
h2l_NKr2_97.5 <- round(boot.ci(results, type ="perc", index = 2)[[4]][1,][5], 6)

res <- data.frame(N, K, P, 
                  NKr2, NKr2_pval, NKr2_2.5, NKr2_97.5,
                  h2l_NKr2, h2l_NKr2_2.5, h2l_NKr2_97.5,
                  tmp[,3:ncol(tmp)])


names(res) <- c("N", "K", "P", 
                "NKr2","NKr2_pval", "NKr2_2.5", "NKr2_97.5", 
                "h2l_NKr2", "h2l_NKr2_2.5", "h2l_NKr2_97.5", 
                "auc1", "auc1_2.5", "auc1_97.5", 
                "auc2", "auc2_2.5", "auc2_97.5","GLM1_Beta","GLM2_Beta")

fwrite(res, file = "AAU_PGC3_SCZ_results.txt", sep = "\t")



##### box plot
# Convert PHENO1 to a factor with labels "Control" and "Case"
prs$PHENO1 <- factor(prs$PHENO1, levels = c(0, 1), labels = c("Control", "Case"))

box.p <- ggplot(prs, aes(x = PHENO1, y = ZSCORE, fill = PHENO1)) +
  geom_boxplot() +
  scale_fill_manual(values = c("#619CFF", "#F8766D")) +
  labs(x = expression(bold("Status")), y = expression(bold("PRS-CS Z Score"))) +
  theme_minimal() +
  theme(
    axis.title = element_text(face = "bold"),
    axis.text.x = element_blank(),  # Remove x-axis labels
    axis.text.y = element_text(size = 14, face = "bold")  # Increase and bold y-axis labels
  )
box.p 

dev.off() 

##### histogram
hist.p <- ggplot(prs, aes(ZSCORE, group = PHENO1, color = PHENO1)) + 
  geom_density(linewidth = 3) + 
  scale_color_manual(values = c("#619CFF", "#F8766D")) + 
  theme_minimal() + 
  theme(
    axis.text.x = element_text(face = "bold", size = 14, color = "black"),
    axis.text.y = element_text(face = "bold", size = 14, color = "black"),
    axis.title = element_text(size = 16, face = "bold")
  ) +
  labs(color = "Status")  # Update legend label

hist.p

dev.off()

