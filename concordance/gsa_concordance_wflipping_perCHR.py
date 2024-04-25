import sys,os
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from scipy.stats import pearsonr
import time


start_time = time.time()

def get_truth_gt(column):
    def calculate_gt(genotype):
        try:
            allele_counts = genotype.split('/')
            if '.' in allele_counts or len(allele_counts) < 2:
                return np.nan
            return sum(int(allele) for allele in allele_counts if allele != '.')
        except ValueError:
            return np.nan

    return column.apply(calculate_gt)


def get_imputed_gt(column):
    def calculate_gt(genotype):
        try:
            allele_counts = genotype.split(':')[0].split("|")
            if '.' in allele_counts or len(allele_counts) < 2:
                return np.nan
            return sum(int(allele) for allele in allele_counts if allele != '.')
        except ValueError:
            return np.nan

    return column.apply(calculate_gt)


def get_imputed_dosages(column):
    def calculate_dosage(genotype):
        try:
            allele_counts = genotype.split(':')[0].split("|")
            if '.' in allele_counts or './.' in genotype:
                return np.nan
            else:
                return float(genotype.split(':')[1])
        except ValueError:
            return np.nan

    return column.apply(calculate_dosage)

def calculate_maf(allele_freq):
    if allele_freq <= 0.5:
        return allele_freq
    else:
        return 1 - allele_freq


def flip_genotypes(genotype_df, reference_df, samples):
    # Merge genotype_df with reference_df on CHROM and POS
    merged_df = pd.merge(genotype_df, reference_df, on=['CHROM', 'POS'], suffixes=('_genotype', '_reference'))
    merged_df = merged_df.set_index('ID', drop=False).rename_axis(None) # setting SNP IDs as row names

    # Identify rows where REF and ALT do not match between the two dataframes
    flipped_rows = merged_df[(merged_df['REF_genotype'] == merged_df['ALT_reference']) & (merged_df['ALT_genotype'] == merged_df['REF_reference'])]

    flipped_mask = genotype_df.index.isin(flipped_rows.index)
    
    # Flip genotypes for mismatched rows
    for sample_id in samples:
        genotype_df.loc[flipped_mask, sample_id] = 2 - genotype_df.loc[flipped_mask, sample_id]
    
    return genotype_df


def compute_metrics(df_imputed, df_true):
    sensitivity = []
    precision = []
    concordance = []

    imputed_array = df_imputed.to_numpy()
    true_array = df_true.to_numpy()

    for i in range(len(imputed_array)):
        alleles_imputed = imputed_array[i]
        alleles_true = true_array[i]
        
        # Find non-NaN indices and common indices
        common_indices = np.intersect1d(np.where(~np.isnan(alleles_imputed)), np.where(~np.isnan(alleles_true)))

        # Filter arrays to common indices
        alleles_imputed = alleles_imputed[common_indices]
        alleles_true = alleles_true[common_indices]

        if np.all(alleles_true == 0):
            sensitivity.append(None)
            precision.append(None)
            concordance.append(None)
            continue

        true_positives = np.sum((alleles_imputed == alleles_true) & (alleles_imputed != 0))
        false_positives = np.sum((alleles_imputed != alleles_true) & (alleles_imputed != 0))
        false_negatives = np.sum((alleles_imputed == 0) & (alleles_true != 0))

        total_positives = true_positives + false_negatives

        sensitivity_snp = true_positives / total_positives if total_positives > 0 else 0
        precision_snp = true_positives / (true_positives + false_positives) if (true_positives + false_positives) > 0 else 0
        concordance_snp = true_positives / (true_positives + false_negatives + false_positives)  if (true_positives + false_negatives + false_positives) > 0 else 0

        sensitivity.append(sensitivity_snp)
        precision.append(precision_snp)
        concordance.append(concordance_snp)

    return sensitivity, precision, concordance


def get_aggregate_R2_generator(geno_df, imp_df, maf_bins=None):
    if maf_bins is None:
        maf_bins = np.linspace(0, 0.5, num=100)

    for i in range(1, len(maf_bins)):
        bin_low, bin_high = maf_bins[i - 1], maf_bins[i]
        maf_cat = f"{bin_low}-{bin_high}"

        idx = np.where((bin_low <= maf_df['MAF']) & (maf_df['MAF'] < bin_high))[0]
        if len(idx) == 0:
            print(f"maf_bin [{bin_low}, {bin_high}] has no SNP! Assigning R2 of NaN")
            yield np.nan, 0, maf_cat
            continue

        truth = np.ravel(geno_df.iloc[idx,])
        imptd = np.ravel(imp_df.iloc[idx,])

        non_missing_idx = np.intersect1d(np.where(~np.isnan(truth))[0], np.where(~np.isnan(imptd))[0])
        my_R2 = pearsonr(truth[non_missing_idx], imptd[non_missing_idx])[0] ** 2
        yield my_R2, len(idx), maf_cat, truth[non_missing_idx], imptd[non_missing_idx]

#get necessary inputs from command line
imputed_filename=sys.argv[1]
gsa_filename=sys.argv[2]
maf_filename=sys.argv[3]
cohort_name=sys.argv[4]

maf_df = pd.read_csv(maf_filename,sep=" ",header=None)
maf_df.columns = ["CHROM","POS","REF","ALT","AF"]
maf_df['SNP'] = (maf_df['CHROM']) + ':' + (maf_df['POS']).astype(str) + ':' + maf_df['REF'].astype(str) + ':' + maf_df['ALT'].astype(str)
maf_df['MAF'] = maf_df['AF'].apply(calculate_maf)
maf_df = maf_df.set_index('SNP', drop=False).rename_axis(None) #setting SNP IDs as row names
print("Successfully loaded and formatted MAF file\n")  

maf_bins = list(np.linspace(0,0.5,num=100))

#read in VCFs as pandas dataframes
df_gsa = pd.read_csv(gsa_filename, sep="\t",comment="#")
df_gsa['ID'] = df_gsa['CHROM'].astype(str) + ':' + df_gsa['POS'].astype(str) + ':' + df_gsa['REF'].astype(str) + ':' + df_gsa['ALT'].astype(str)
df_gsa = df_gsa.set_index('ID', drop=False).rename_axis(None) #setting SNP IDs as row names
print("Successfully loaded GSA VCF file\n")

df_imputed = pd.read_csv(imputed_filename, sep="\t",comment="#")
df_imputed['ID'] = (df_imputed['CHROM']) + ':' + (df_imputed['POS']).astype(str)  + ':' + df_imputed['REF'].astype(str) + ':' + df_imputed['ALT'].astype(str)
df_imputed = df_imputed.set_index('ID', drop=False).rename_axis(None)
print("Successfully loaded imputed VCF file\n")


#initially identify SNPs by CHR:POS:REF:ALT and CHR:POS:ALT:REF to collect all possible alleles
gsa_snps = df_gsa['CHROM'].astype(str) + ':' + df_gsa['POS'].astype(str) + ':' + df_gsa['REF'].astype(str) + ':' + df_gsa['ALT'].astype(str)
gsa_snps_flip =  df_gsa['CHROM'].astype(str) + ':' + df_gsa['POS'].astype(str) + ':' + df_gsa['ALT'].astype(str) + ':' + df_gsa['REF'].astype(str)
gsa_snps = pd.concat([gsa_snps, gsa_snps_flip], axis=0)


imputed_snps = (df_imputed['CHROM']) + ':' + (df_imputed['POS']).astype(str)  + ':' + df_imputed['REF'].astype(str) + ':' + df_imputed['ALT'].astype(str)
imputed_snps_flip = (df_imputed['CHROM']) + ':' + (df_imputed['POS']).astype(str)  + ':' + df_imputed['ALT'].astype(str) + ':' + df_imputed['REF'].astype(str)
imputed_snps = pd.concat([imputed_snps, imputed_snps_flip],axis=0)

#Combine SNP lists to get all SNPs with an exact or flipped match (this is where non-matched multiallelic SNPs get dropped)
initial_snps = set(imputed_snps).intersection(set(gsa_snps))
imputed_snps_indices = df_imputed.index.intersection(initial_snps)
gsa_snps_indices = df_gsa.index.intersection(initial_snps)

#Get intersecting samples between the two VCFs
intersecting_columns = df_gsa.columns.intersection(df_imputed.columns)
samples = intersecting_columns.drop(['CHROM', 'POS', 'ID', 'REF', 'ALT', 'QUAL', 'FILTER', 'INFO', 'FORMAT'])

#get initial subset of dataframes                                                           
imp_intersect = df_imputed.loc[imputed_snps_indices, intersecting_columns]
gsa_intersect = df_gsa.loc[gsa_snps_indices, intersecting_columns]
maf_df = maf_df.loc[maf_df['SNP'].isin(initial_snps)]

#########
n_imptd_snps = len(imputed_snps_indices)
n_gsa_snps = len(gsa_snps_indices)
n_samples = len(samples)

print(f"Subsetting of VCF files includes:\n {n_gsa_snps} GSA SNPs\n {n_imptd_snps} Imputed SNPs\n {n_samples} Samples\n")
########


gsa_matrix = gsa_intersect.apply(lambda col: get_truth_gt(col) if col.name in samples else col)
print("Successfully extracted GSA genotypes from VCF\n")
        
imputed_ds_matrix= imp_intersect.apply(lambda col: get_imputed_dosages(col) if col.name in samples else col)
print("Successfully extracted imputed dosages from VCF\n")

imputed_gt_matrix = imp_intersect.apply(lambda col: get_imputed_gt(col) if col.name in samples else col)
print("Successfully extracted imputed genotypes from VCF\n")


#resetting SNP IDs to CHR:POS to aid in identifying flipping genotypes 
gsa_matrix['ID'] = gsa_matrix['CHROM'].astype(str) + ':' + gsa_matrix['POS'].astype(str)
gsa_matrix = gsa_matrix.set_index('ID', drop=False).rename_axis(None) #setting SNP IDs as row names

imputed_ds_matrix['ID'] = imputed_ds_matrix['CHROM'].astype(str) + ':' + imputed_ds_matrix['POS'].astype(str)
imputed_ds_matrix = imputed_ds_matrix.set_index('ID', drop=False).rename_axis(None) #setting SNP IDs as row names

imputed_gt_matrix['ID'] = imputed_gt_matrix['CHROM'].astype(str) + ':' + imputed_gt_matrix['POS'].astype(str)
imputed_gt_matrix = imputed_gt_matrix.set_index('ID', drop=False).rename_axis(None) #setting SNP IDs as row names

maf_df['SNP'] = (maf_df['CHROM']) + ':' + (maf_df['POS']).astype(str)
maf_df = maf_df.set_index('SNP', drop=False).rename_axis(None) #setting SNP IDs as row names


##### flip genotypes to match allele orientation between the VCFs and MAF dataframe
imputed_gt_fixed = flip_genotypes(imputed_gt_matrix, maf_df, samples)
print("Succesfully flipped alleles in imputed genotype dataframe\n")

imputed_ds_fixed = flip_genotypes(imputed_ds_matrix, maf_df, samples)
print("Succesfully flipped alleles in imputed dosage dataframe\n")

#note if MAF file is derived from GSA allele frequencies, this step can be skipped
gsa_fixed = flip_genotypes(gsa_matrix, maf_df, samples)
print("Successfully flipped alleles in GSA dataframe\n") 


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


#redo SNP intersection after allele flipping                    
intersecting_snps = gsa_fixed.index.intersection(imputed_gt_fixed.index)
gsa_fixed = gsa_fixed.loc[intersecting_snps, samples]
imputed_ds_fixed = imputed_ds_fixed.loc[intersecting_snps, samples]
imputed_gt_fixed = imputed_gt_fixed.loc[intersecting_snps, samples]
maf_df = maf_df.loc[maf_df['SNP'].isin(intersecting_snps)]

#check to make sure number of SNPs and number of samples match
if len(imputed_gt_fixed.index) != len(gsa_fixed.index):
    sys.exit("The number of SNPS in the truth and imputated datasets are not equal\n")
if len(imputed_gt_fixed.columns) != len(gsa_fixed.columns):
    sys.exit("The number of samples in the truth and imputated datasets are not equal\n")

print("Number of samples and SNPs now equal between GSA and imputed dataframes\n")

#make sure samples are in the same order
gsa_fixed = gsa_fixed.reindex(columns=samples)
imputed_gt_fixed = imputed_gt_fixed.reindex(columns=samples)
imputed_ds_fixed = imputed_ds_fixed.reindex(columns=samples)


########## compute non-ref concordance ############
print("Computing non-reference concordance ... \n")
s, p, c = compute_metrics(imputed_gt_fixed, gsa_fixed)

print("Computed non-reference concordance per SNP\n")
df_metrics = pd.DataFrame({'SNP': imputed_gt_fixed.index, 'Sensitivity': s, 'Precision': p, 'Non-Ref Concordance': c})
df_metrics = df_metrics.merge(maf_df,on=["SNP"])
df_metrics_fname = cohort_name + "_nonref_concordance_table.csv"
df_metrics.to_csv(df_metrics_fname,index=False)
print(f"Per-SNP non-reference concordance CSV written to {df_metrics_fname}\n")

df_metrics['MAF_bin'] = pd.cut(df_metrics['MAF'], maf_bins)
df_filtered = df_metrics[df_metrics['Non-Ref Concordance'].notna()]
df_binned = df_filtered.groupby('MAF_bin').agg({
    'Non-Ref Concordance': 'mean',
    'SNP': 'count'}).reset_index()

df_binned.columns = ["MAF_bin","Average Non-Ref Concordance", "Number of SNPs"]

df_binned_fname = cohort_name + "_binned_nonref_concordance_table.csv"
df_binned.to_csv(df_binned_fname,index=False)
print(f"Binned non-reference concordance CSV written to {df_binned_fname}\n")


########### compute aggregate R2 ############
print("Computing aggregated R2...\n")
r2_list, nsnps_list, maf_list, true_list, imptd_list = zip(*get_aggregate_R2_generator(gsa_fixed, imputed_ds_fixed))
print("Computed aggregated R2\n")
aggregation_df = pd.DataFrame({'maf_category': maf_list, 'number_of_snps': nsnps_list, 'R2': r2_list})
gt_lists = pd.DataFrame(list(zip(maf_list,true_list, imptd_list)), columns=['MAF', 'True_list','Imputed_list'])
gt_list_name = cohort_name + "_genotypes_per_MAFbin_table.tsv"
gt_lists.to_csv(gt_list_name,index=False,sep="\t")


agg_df_name = cohort_name + "_aggregation_R2_table.csv"
aggregation_df.to_csv(agg_df_name,index=False)
print(f"Aggregated R2 CSV written to {agg_df_name}\n")

end_time = time.time()  # Record the end time
execution_time = end_time - start_time  # Calculate the execution time#
print(f"Execution time: {execution_time} seconds")
print("Script completed! Exiting.\n")

