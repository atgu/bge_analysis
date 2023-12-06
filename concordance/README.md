# Evaluating imputation accuracy

We used a few metrics for evaluating GLIMPSE's imputation accuracy:
1. Concordance (including sensitivity, precision, and non-reference concordance)
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

## Sensitivity, precision, and non-reference concordance

```shell
julia concordance.jl --truth file1 --impt file2 --out file3
```
Required arguments:
+ `--truth`: Full file path to the ground truth genotypes (must end in `.vcf` or `.vcf.gz`), or a file containing a list of VCF files
+ `--impt`: Full file path to the imputed genotypes (must end in `.vcf` or `.vcf.gz`), or a file containing a list of VCF files
+ `--out`: Full file path to output file name. Output will be non-reference concordance (i.e the percent of non-reference alleles that were correctly imputed)

Optional arguments include:
+ `--summary`: An optional which allows users to NOT consider certain SNPs that exist in ground truth or imputed genotype files. Must be either a full file path to a comma or tab separated summary file (header must include at least 5 columns with the names CHR/POS/REF/ALT/isImputed), or a file containing a list of summary files. The scripts are set up so that only SNPs listed as `true` in the `isImputed` column will be considered. 
+ `--maf-bins`: Comma-separated list of minor allele frequencies used to bin SNPs. Defaults to `0.0,0.0005,0.001,0.004,0.0075,0.0125,0.04,0.1,0.2,0.5`

### Details of sensitivity, precision, and non-reference concordances

For each SNP, we can create the following contingency table:

|       |         | Imputed            | Imputed            |
|-------|---------|--------------------|--------------------|
|       |         | ALT (1)            | REF (0)            |
| truth | ALT (1) | true positive (A)  | false negative (B) |
| truth | REF (0) | false positive (C) | false negative (D) |

Let 1 = ALT and 0 = REF, then we have the following cases

+ truth = 0 and imputed = 0: D += 2
+ truth = 1 and imputed = 0: D += 1 and B += 1
+ truth = 2 and imputed = 0: B += 1
+ truth = 0 and imputed = 1: D += 1 and C += 1
+ truth = 1 and imputed = 1: A += 1 and D += 1
+ truth = 2 and imputed = 1: A += 1 and B += 1
+ truth = 0 and imputed = 2: C += 2
+ truth = 1 and imputed = 2: A += 1 and C += 1
+ truth = 2 and imputed = 2: A += 2

Then

+ `Sensitivity` = true positives / (true positives + false negaitves) = A / (A+B). (NOTE that the total number of ALT alleles in the array data is A+B)
+ `Precision` = true positives / (true positives + false positives) = A / (A+C)
+ `Non-reference concordance` = A / (A+B+C). To see why, consider the denominator to be the total number of ALT alleles called the *either* array/truth dataset or imputed dataset, i.e. (A+B+C). The numerator is where they agree, i.e. A.

## Aggregate R2

```shell
julia aggregate_r2.jl --truth file1 --impt file2 --out file3
```
Required arguments:
+ `--truth`: Full file path to the ground truth genotypes (must end in .vcf or .vcf.gz), or a file containing a list of VCF files
+ `--impt`: Full file path to the imputed genotypes (must end in .vcf or .vcf.gz), or a file containing a list of VCF files
+ `--out`: Full file path to output file name. Output will be non-reference concordance (i.e the percent of non-reference alleles that were correctly imputed)

Optional arguments include:
+ `--summary`: An optional which allows users to NOT consider certain SNPs that exist in ground truth or imputed genotype files. Must be either a full file path to a comma or tab separated summary file (header must include at least 5 columns with the names CHR/POS/REF/ALT/isImputed), or a file containing a list of summary files. The scripts are set up so that only SNPs listed as `true` in the `isImputed` column will be considered. 
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
+ `--summary`: An optional which allows users to NOT consider certain SNPs that exist in ground truth or imputed genotype files. Must be either a full file path to a comma or tab separated summary file (header must include at least 5 columns with the names CHR/POS/REF/ALT/isImputed), or a file containing a list of summary files. The scripts are set up so that only SNPs listed as `true` in the `isImputed` column will be considered. 
+ `--maf-bins`: Comma-separated list of minor allele frequencies used to bin SNPs. Defaults to `0.0,0.0005,0.001,0.004,0.0075,0.0125,0.04,0.1,0.2,0.5`

## Plotting

+ Plotting code are provided in `plots.jl`. 
+ See `examples.ipynb` for plotting examples. 

## Memory Requirements and runtime

These scripts import the imputed and ground truth genotypes into double-precision matrices (without multithreading). Thus, we *highly* recommend one to separate data by chromosome, and filter the imputed data so that they are on roughly the same set of SNPs as the ground truth data. 

As a rule of thumb, the provided scripts should take at most an hour on ~1000 samples and ~400k SNPs. A progress bar will be automatically displayed for routines such as importing data. 

## Bugs, usage difficulties, and feature requests?

Please submit an issue for any problems and I will take a look asap. Alternatively, feel free to reach out to Benjamin Chu (bbchu@stanford.edu). 
