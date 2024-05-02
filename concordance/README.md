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
