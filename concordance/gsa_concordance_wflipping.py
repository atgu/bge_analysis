import sys,os
import pandas as pd
import numpy as np 
import matplotlib.pyplot as plt
from scipy.stats import pearsonr


########## functions ###########

def get_truth_gt(genotypes):
    allele_counts = genotypes.str.split('/', expand=True)
    gt = allele_counts.apply(lambda x: np.nan if '.' in x.values else x.astype(int).sum(), axis=1)
    return gt


def get_imputed_gt(genotypes):
    allele_counts = genotypes.str.split(':').str[0].str.split("|")
    gt = allele_counts.apply(lambda x: int(x[0]) + int(x[1]))
    return gt


def get_imputed_dosages(genotypes):
    # Extract the allele counts from the imputed-style genotype string
    dosages = genotypes.str.split(':').str[1]
    dosages = dosages.astype(float)
    return dosages


def calculate_maf(allele_freq):
    if allele_freq <= 0.5:
        return allele_freq
    else:
        return 1 - allele_freq


def flip_genotypes(genotype_df, reference_df, samples):
    # Merge genotype_df with reference_df on CHROM and POS
    merged_df = pd.merge(genotype_df, reference_df, on=['CHROM', 'POS'], suffixes=('_genotype', '_reference'))

    # Identify rows where REF and ALT do not match between the two dataframes
    mismatch_rows = merged_df[(merged_df['REF_genotype'] != merged_df['REF_reference']) | (merged_df['ALT_genotype'] != merged_df['ALT_reference'])]

    # Flip genotypes for mismatched rows
    for sample_id in samples:
        genotype_df[sample_id] = genotype_df.apply(lambda row: 2 - row[sample_id] if row.name in mismatch_rows.index else row[sample_id], axis=1)

    # Reset index before dropping rows
    genotype_df.reset_index(drop=True, inplace=True)

    # Remove rows where the mismatch still exists after flipping
    genotype_df = genotype_df.drop(mismatch_rows.index)

    return genotype_df


def compute_metrics(df_imputed, df_true):
    sensitivity = []
    precision = []
    concordance = []
    
    for idx in df_imputed.index:
        alleles_imputed = df_imputed.loc[idx]
        alleles_true = df_true.loc[idx]
        
        true_positives = sum((alleles_imputed == alleles_true) & (alleles_imputed != 0))  # Count non-reference alleles that are correctly imputed
        false_positives = sum((alleles_imputed != alleles_true) & (alleles_imputed != 0))  # Count non-reference alleles incorrectly imputed
        false_negatives = sum((alleles_imputed == 0) & (alleles_true != 0))  # Count non-reference alleles missed by imputation
        
        total_positives = true_positives + false_negatives

        sensitivity_snp = true_positives / total_positives if total_positives > 0 else 0
        precision_snp = true_positives / (true_positives + false_positives) if (true_positives + false_positives) > 0 else 0
        
        concordance_snp = true_positives / (true_positives + false_negatives + false_positives)  if (true_positives + false_negatives + false_positives) > 0 else 0

        sensitivity.append(sensitivity_snp)
        precision.append(precision_snp)
        concordance.append(concordance_snp)
    
    return sensitivity, precision, concordance


def get_aggregate_R2(geno_df, imp_df, maf_bins=None):

    if maf_bins is None:
        maf_bins = list(np.linspace(0,0.5,num=100))

    R2 = []
    nsnps = []
    maf_category = []

    for i in range(1, len(maf_bins)):
        # find SNPs in current bin
        bin_low, bin_high = maf_bins[i-1], maf_bins[i]
        idx = np.where((bin_low <= maf_df['MAF']) & (maf_df['MAF'] < bin_high))[0]
        maf_cat = str(bin_low) + "-" + str(bin_high)
        maf_category.append(maf_cat)

        if len(idx) == 0:
            print(f"maf_bin [{bin_low}, {bin_high}] has no SNP! Assigning R2 of NaN")
            R2.append(np.nan)
            nsnps.append(0)
            continue
        
        # compute squared pearson correlation
        truth = np.ravel(geno_df.iloc[idx,])
        imptd = np.ravel(imp_df.iloc[idx,])

        non_missing_idx = np.intersect1d(np.where(~np.isnan(truth))[0], np.where(~np.isnan(imptd))[0])

        my_R2 = pearsonr(truth[non_missing_idx], imptd[non_missing_idx])[0] ** 2
        R2.append(my_R2)
        nsnps.append(len(idx))

    return R2, nsnps, maf_category


#######################################################################################
        
#get necessary inputs from command line
imputed_filename=sys.argv[1]
gsa_filename=sys.argv[2]
maf_filename=sys.argv[3]
cohort_name=sys.argv[4]

# read in allele frequency info (tab separated file with columns CHROM POS REF ALT AF
maf_df = pd.read_csv(maf_filename,sep=" ",header=None)
maf_df.columns = ["CHROM","POS","REF","ALT","AF"]
maf_df['SNP'] = (maf_df['CHROM']) + ':' + (maf_df['POS']).astype(str)
maf_df['MAF'] = maf_df['AF'].apply(calculate_maf)
maf_df = maf_df.set_index('SNP', drop=False).rename_axis(None) #setting SNP IDs as row names

maf_bins = list(np.linspace(0,0.5,num=100))
maf_bins2 = [0.0, 0.0005, 0.001, 0.005, 0.01, 0.02, 0.03, 0.04, 0.05, 0.1, 0.2, 0.3, 0.4, 0.5]

#read in as pandas dataframes
df_imputed = pd.read_csv(imputed_filename, sep="\t",comment="#")
df_gsa = pd.read_csv(gsa_filename, sep="\t",comment="#")


#format dataframes
df_gsa['ID'] = df_gsa['CHROM'].astype(str) + ':' + df_gsa['POS'].astype(str)
df_gsa = df_gsa.set_index('ID', drop=False).rename_axis(None) #setting SNP IDs as row names

df_imputed['ID'] = (df_imputed['CHROM']) + ':' + (df_imputed['POS']).astype(str)
df_imputed = df_imputed.set_index('ID', drop=False).rename_axis(None)  


intersecting_snps = df_gsa.index.intersection(df_imputed.index)
intersecting_columns = df_gsa.columns.intersection(df_imputed.columns)
samples = intersecting_columns.drop(['CHROM', 'POS', 'ID', 'REF', 'ALT', 'QUAL', 'FILTER', 'INFO', 'FORMAT'])

#get initial subset of dataframes 
imp_intersect = df_imputed.loc[intersecting_snps, intersecting_columns]
gsa_intersect = df_gsa.loc[intersecting_snps, intersecting_columns]
maf_df = maf_df.loc[maf_df['SNP'].isin(intersecting_snps)]

    
# get genotype information in correct format
imputed_gt_matrix = imp_intersect.apply(lambda col: get_imputed_gt(col) if col.name in samples else col)
imputed_ds_matrix = imp_intersect.apply(lambda col: get_imputed_dosages(col)  if col.name in samples else col)
gsa_matrix = gsa_intersect.apply(lambda col: get_truth_gt(col) if col.name in samples else col)


#fix the dataframes by flipping alleles as needed (and remove SNPs that still have non-matching REF/ALT alleles even after flipping)
imputed_gt_fixed = flip_genotypes(imputed_gt_matrix, maf_df, samples)
imputed_ds_fixed = flip_genotypes(imputed_ds_matrix, maf_df, samples)
gsa_fixed = flip_genotypes(gsa_matrix, maf_df, samples)


#resetting indices of fixed dataframes
imputed_gt_fixed['ID'] = imputed_gt_fixed['CHROM'].astype(str) + ':' + imputed_gt_fixed['POS'].astype(str)
imputed_gt_fixed = imputed_gt_fixed.set_index('ID', drop=False).rename_axis(None) #setting SNP IDs as row names
imputed_gt_fixed = imputed_gt_fixed.iloc[:, 9 :] #remove snp data from df

imputed_ds_fixed['ID'] = imputed_ds_fixed['CHROM'].astype(str) + ':' + imputed_ds_fixed['POS'].astype(str)
imputed_ds_fixed = imputed_ds_fixed.set_index('ID', drop=False).rename_axis(None) #setting SNP IDs as row names
imputed_ds_fixed = imputed_ds_fixed.iloc[:, 9 :] #remove snp data from df

gsa_fixed['ID'] = gsa_fixed['CHROM'].astype(str) + ':' + gsa_fixed['POS'].astype(str)
gsa_fixed = gsa_fixed.set_index('ID', drop=False).rename_axis(None) #setting SNP IDs as row names
gsa_fixed = gsa_fixed.iloc[:, 9 :] #remove snp data from df


#redo SNP intersection after allele flipping / removing rows 
intersecting_snps = gsa_fixed.index.intersection(imputed_gt_fixed.index)
gsa_fixed = gsa_fixed.loc[intersecting_snps, samples]
imputed_ds_fixed = imputed_ds_fixed.loc[intersecting_snps, samples]
imputed_gt_fixed = imputed_gt_fixed.loc[intersecting_snps, samples]
maf_df = maf_df.loc[maf_df['SNP'].isin(intersecting_snps)]

#check to make sure number of SNPs and number of samples match
if len(imputed_gt_fixed.index) != len(gsa_fixed.index):
    print("The number of SNPS in the truth and imputated datasets are not equal")
    #break out of script
    
if len(imputed_gt_fixed.columns) != len(gsa_fixed.columns):
    print("The number of samples in the truth and imputated datasets are not equal")
    #break out of script                                                                                                                             

#make sure samples are in the same order
gsa_fixed = gsa_fixed.reindex(columns=samples)
imputed_gt_fixed = imputed_gt_fixed.reindex(columns=samples)
imputed_ds_fixed = imputed_ds_fixed.reindex(columns=samples)



########## compute non-ref concordance ############
s, p, c = compute_metrics(imputed_gt_fixed, gsa_fixed)
df_metrics = pd.DataFrame({'SNP': imputed_gt_fixed.index, 'Sensitivity': s, 'Precision': p, 'Non-Ref Concordance': c})


df_metrics = df_metrics.merge(maf_df,on=["SNP"])
df_metrics_fname = cohort_name + "_nonref_concordance_table.csv"
df_metrics.to_csv(df_metrics_fname,index=False)


df_metrics['MAF_bin'] = pd.cut(df_metrics['MAF'], maf_bins)

df_binned = df_metrics.groupby('MAF_bin').agg({
    'Non-Ref Concordance': 'mean',
    'SNP': 'count'
}).reset_index()

df_binned.columns = ["MAF_bin","Average Non-Ref Concordance", "Number of SNPs"]

df_binned_fname = cohort_name + "_binned_nonref_concordance_table.csv"
df_binned.to_csv(df_binned_fname,index=False)

###### for emphasizing rarer SNPs ########
df_metrics['MAF_bin2'] = pd.cut(df_metrics['MAF'], maf_bins2)

df_binned2 = df_metrics.groupby('MAF_bin2').agg({
    'Non-Ref Concordance': 'mean',
    'SNP': 'count'
}).reset_index()

df_binned2.columns = ["MAF_bin","Average Non-Ref Concordance", "Number of SNPs"]

df_binned2_fname = cohort_name + "_rare_binned_nonref_concordance_table.csv"
df_binned2.to_csv(df_binned2_fname,index=False)



#### plotting non-ref concordance
#plt.scatter(df_binned['MAF_bin'].astype(str), df_binned['Average Non-Ref Concordance'],df_binned['Number of SNPs'].astype(float)/100)
#plt.xlabel('MAF Bin')
#plt.ylabel('Average Non-Ref Concordance')
#plt.title('Average Non-Ref Concordance by MAF Bin')

# Set specific tick marks for MAF bins
#tick_positions = range(0, len(df_binned['MAF_bin']), len(df_binned['MAF_bin']) // 5)
#plt.xticks(ticks=tick_positions, labels=['0.0', '0.1', '0.2', '0.3', '0.4', '0.5'])
#plt.ylim(0.0, 1.0)
#non_ref_fig_name = cohort_name + "_non-reference_concordance_figure.pdf"
#plt.savefig(non_ref_fig_name, format="pdf")
#plt.show()


########### compute aggregate R2 ############
r2, nsnps, maf = get_aggregate_R2(gsa_fixed,imputed_ds_fixed)
aggregation_df = pd.DataFrame(list(zip(maf, nsnps, r2)), columns=['maf_category', 'number_of_snps', 'R2'])

agg_df_name = cohort_name + "_aggregation_R2_table.csv"
aggregation_df.to_csv(agg_df_name,index=False)



########## for emphasizing rarer snps #######
r2_2, nsnps_2, maf_2 = get_aggregate_R2(gsa_fixed,imputed_ds_fixed,maf_bins2)
aggregation_df2 = pd.DataFrame(list(zip(maf_2, nsnps_2, r2_2)), columns=['maf_category', 'number_of_snps', 'R2'])

agg_df2_name = cohort_name + "_rare_aggregation_R2_table.csv"
aggregation_df2.to_csv(agg_df2_name,index=False)



#### plotting aggregate R2
#plt.scatter(aggregation_df['maf_category'], aggregation_df['R2'],s=aggregation_df['number_of_snps'].astype(float)/100)
#plt.xlabel('MAF Bin')
#plt.ylabel('Aggregate R2')
#plt.title('Aggregate R2 by MAF Bin')

# Set specific tick marks for MAF bins
#tick_positions = range(0, len(aggregation_df['maf_category']), len(aggregation_df['maf_category']) // 5)
#plt.xticks(ticks=tick_positions, labels=['0.0', '0.1', '0.2', '0.3', '0.4', '0.5'])
#plt.ylim(0.0, 1.0)

#agg_fig_name = cohort_name + "_aggregation_R2_figure.pdf"
#plt.savefig(agg_fig_name, format="pdf")
#plt.show()
