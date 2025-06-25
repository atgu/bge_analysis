We used two main metrics for evaluating GLIMPSE's imputation accuracy:

1. Concordance (including sensitivity, precision, and non-reference concordance)
2. Aggregate R2 (squared Pearson's correlation for all SNPs in a given MAF bin)

How is concordance calculated? 

For each SNP, we can create the following contingency table:

|            |              | **Imputed**        | **Imputed**        |
|------------|--------------|--------------------|--------------------|
|            |              | ALT (1)            | REF (0)            |
| **truth**  | ALT (1)      | true positive (A)  | false negative (B) |
| **truth**  | REF (0)      | false positive (C) | true negative (D)  |

Let 1 = ALT and 0 = REF, then we have the following cases:

truth = 0 and imputed = 0: D += 2  
truth = 1 and imputed = 0: D += 1 and B += 1  
truth = 2 and imputed = 0: B += 1  
truth = 0 and imputed = 1: D += 1 and C += 1  
truth = 1 and imputed = 1: A += 1 and D += 1  
truth = 2 and imputed = 1: A += 1 and B += 1  
truth = 0 and imputed = 2: C += 2  
truth = 1 and imputed = 2: A += 1 and C += 1  
truth = 2 and imputed = 2: A += 2  

Then

`Sensitivity` = true positives / (true positives + false negatives) = A / (A+B). (NOTE that the total number of ALT alleles in the array data is A+B)  
`Precision` = true positives / (true positives + false positives) = A / (A+C)  
`Non-reference concordance` = A / (A+B+C). To see why, consider the denominator to be the total number of ALT alleles called the either array/truth dataset or imputed dataset, i.e. (A+B+C). The numerator is where they agree, i.e. A.  

This script requires the following python packages to be installed:
1. Pandas
2. Numpy
3. Matplotlib
4. Scipy.stats

Four command-line arguments are required: 
1. Imputed VCF file (can be .vcf.gz or .vcf)
2. QC’ed GSA VCF file (can be .vcf.gz or .vcf)
3. MAF file
4. Cohort name (used for output file prefix) 


Run the script via:

`python gsa_concordance.py <imputed.vcf> <gsa.vcf> <maf.txt> <cohort>`

*Note for VCF files, in order for pandas to read in the column headers for the dataframe, the “#CHROM” needs to be edited to just “CHROM”, since any lines that begin with “#” are considered comments and ignored. The following are the bcftools commands we ran to format the imputed NeuroGAP data: 

#first, update the sample names to match those in the GSA VCF, and index                                                                                 
`bcftools reheader -s {sample_list} -o {name}_temp.bcf {input_vcf}`                                                                                                               
`bcftools index {name}_temp.bcf`

#next keep only variants that are in the GSA array                                                                                            
`bcftools view -R {variant_list} -Oz -o {name}_temp2.vcf.gz {name}_temp.bcf`

#unzip                                                                                                                                        
`gunzip {name}_temp2.vcf.gz`

#remove “#” from CHROM header line so it can be loaded as pandas dataframe                                                                      
`sed 's/#CHROM/CHROM/g' {name}_temp2.vcf > {name}_temp3.vcf`

#gzip and write to a standalone "vcf" file without indexing (won't work with bcftools because of the missing #CHROM)                                           
`gzip -c {name}temp3.vcf > {output_vcf}`


The MAF file should be space-delimited with 5 columns: CHR, POS, REF, ALT, and AF (we convert to MAF from AF in the python script itself). We generated this file with the following bcftools commands:

`bcftools norm -d all -o dedup.bcf <GSA.vcf>`

`bcftools index dedup.bcf`

`bcftools +fill-tags dedup.bcf -Ob -o out.bcf`

`bcftools query -f '%CHROM %POS %REF %ALT %AF\n' out.bcf > <GSA_allele_freqs.txt>`


For ~150 individuals and ~400k SNPs, the script takes about ~3min to run locally, with 16GB RAM across 4 cores. 
For ~1,000 individuals and ~400k SNPs, the script takes about 1 hour to run, with 64GB of RAM across 4 cores. 


For concordance and aggregated R2 computed per local ancestry haplotype, please see the scripts in the `julia` directory. 
