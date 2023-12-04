# Evaluating imputation accuracy

We used a few metrics for evaluating GLIMPSE's imputation accuracy:
1. Non-reference concordance
2. Aggregate R2 (squared Pearson's correlation)
3. Aggregate R2 but stratefying SNPs by their local ancestry backgrounds. 

See `examples.ipynb` for how to run these scripts and expected outputs. 

## Installation

To run the scripts, 

1. Download and install [Julia](https://julialang.org/downloads/) (any version above 1.6.7 should work)
2. Within Julia, install required packages by 
```julia
using Pkg
Pkg.add(["VCFTools", "DataFrames", "CSV", "DelimitedFiles", "Statistics", "ProgressMeter", "ArgParse", "Plots"])
```
3. Obtain a copy of this repository via
```
git clone https://github.com/biona001/bge_analysis.git
```
4. A few scripts is located in `bge_analysis/concordance/scripts`. Follow the examples below to run them for each analysis

## Non-reference concordance

```shell
julia concordance.jl --truth file1 --impt file2 --out file3
```
Required arguments:
+ `--truth`: Full file path to the ground truth genotypes (must end in `.vcf` or `.vcf.gz`), or a file containing a list of VCF files
+ `--impt`: Full file path to the imputed genotypes (must end in `.vcf` or `.vcf.gz`), or a file containing a list of VCF files
+ `--out`: Full file path to output file name. Output will be non-reference concordance (i.e the percent of non-reference alleles that were correctly imputed)

Optional arguments include:
+ `--summary`: An optional input. Either Full file path to a comma or tab separated summary file (header must include at least 6 columns with the names CHR/POS/REF/ALT/AF/isImputed), or a file containing a list of summary files
+ `--maf-bins`: Comma-separated list of minor allele frequencies used to bin SNPs. Defaults to `0.0,0.0005,0.001,0.004,0.0075,0.0125,0.04,0.1,0.2,0.5`

## Aggregate R2

```shell
julia aggregate_r2.jl --truth file1 --impt file2 --out file3
```
Required arguments:
+ `--truth`: Full file path to the ground truth genotypes (must end in .vcf or .vcf.gz), or a file containing a list of VCF files
+ `--impt`: Full file path to the imputed genotypes (must end in .vcf or .vcf.gz), or a file containing a list of VCF files
+ `--out`: Full file path to output file name. Output will be non-reference concordance (i.e the percent of non-reference alleles that were correctly imputed)

Optional arguments include:
+ `--summary`: An optional input. Either Full file path to a comma or tab separated summary file (header must include at least 6 columns with the names CHR/POS/REF/ALT/AF/isImputed), or a file containing a list of summary files
+ `--maf-bins`: Comma-separated list of minor allele frequencies used to bin SNPs. Defaults to `0.0,0.0005,0.001,0.004,0.0075,0.0125,0.04,0.1,0.2,0.5`

## Aggregate R2 based for different local ancestry backgrounds

Note that
+ Before computing local ancestries, a standard step is to phase the target dataset. We recommend first running the [conform_gt](https://faculty.washington.edu/browning/conform-gt.html) program and then phase with [Beagle 5.4](https://faculty.washington.edu/browning/beagle/beagle.html) 
+ We assume that local ancestries were computed on the ground truth data via the [Rfmix2](https://github.com/slowkoni/rfmix/blob/master/MANUAL.md) software, outputting a `.msp.tsv` file. 

```shell
julia local_ancestry_r2.jl --truth file1 --impt file2 --msp file3 --out file4
```

Required arguments:
+ `--truth`: Full file path to the ground truth genotypes (must end in `.vcf` or `.vcf.gz`), or a file containing a list of VCF files
+ `--impt`: Full file path to the imputed genotypes (must end in `.vcf` or `.vcf.gz`), or a file containing a list of VCF files
+ `--msp`: Full file path to the `.msp.tsv` file (output of rfmix2), or a file containing a list of `.msp.tsv` files
+ `--out`: Full file path to output file name. Output will be non-reference concordance (i.e the percent of non-reference alleles that were correctly imputed)

Optional arguments include:
+ `--summary`: An optional input. Either Full file path to a comma or tab separated summary file (header must include at least 6 columns with the names CHR/POS/REF/ALT/AF/isImputed), or a file containing a list of summary files
+ `--maf-bins`: Comma-separated list of minor allele frequencies used to bin SNPs. Defaults to `0.0,0.0005,0.001,0.004,0.0075,0.0125,0.04,0.1,0.2,0.5`

## Plotting

+ Plotting code are provided in `plots.jl`. 
+ See `examples.ipynb` for plotting examples. 

## Memory Requirements and runtime

These scripts import the imputed and ground truth genotypes into numeric matrices (without multithreading). Thus, we *highly* recommend one to separate data by chromosome, and filter the imputed data so that they are on roughly the same set of SNPs as the ground truth data. 


## Bugs, usage difficulties, and feature requests

Please submit an issue for any problems and I will take a look asap. Alternatively, feel free to reach out to Benjamin Chu (bbchu@stanford.edu). 
