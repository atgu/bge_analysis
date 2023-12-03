# usage: julia local_ancestry_r2.jl input1 input2 input3 input4
# input1: Full file path to the ground truth genotypes (must end in .vcf or .vcf.gz), 
#         or a file containing a list of VCF files
# input2: Full file path to the imputed genotypes (must end in .vcf or .vcf.gz), 
#         or a file containing a list of VCF files
# input3: Full file path to a comma or tab separated summary file (header must 
#         include at least 6 columns with the names CHR/POS/REF/ALT/AF/isImputed), 
#         or a file containing a list of summary files
# input4: Full file path to the .msp.tsv file (output of rfmix2), or a file 
#         containing a list of .msp.tsv files
# input5: Number of ancestries used (should be an integer)
# input6: Full file path to output file name. Output will be aggregate R2 for 
#         each ancestry for different MAF bins. If a list of files were given in 
#         inputs 1-3, then we will average the aggregate R2 for each file (but
#         accounting for the number of SNPs in each file)

include("utilities.jl")

# inputs to this script
input1 = ARGS[1]
input2 = ARGS[2]
input3 = ARGS[3]
input4 = ARGS[4]
input5 = parse(Int, ARGS[5])
outfile = ARGS[6]

# default minor allele frequency bins and ancestry names
maf_bins = [0.0, 0.0005, 0.00100, 0.004, 0.0075, 0.0125, 0.04, 0.1, 0.2, 0.5]
ancestry_names = ["AFR", "AMR", "EUR"]

# compute non-reference concordance
@info "Computing aggregate R2 by local ancestries"
input1_isvcf = endswith(input1, ".vcf") || endswith(input1, ".vcf.gz")
input2_isvcf = endswith(input2, ".vcf") || endswith(input2, ".vcf.gz")
if input1_isvcf && input2_isvcf
    genotype_files = input1
    imputed_files = input2
    summary_files = input3
    msp_files = input4
elseif !input1_isvcf && !input2_isvcf
    genotype_files = readdlm(input1) |> vec |> Vector{String}
    imputed_files = readdlm(input2) |> vec |> Vector{String}
    summary_files = readdlm(input3) |> vec |> Vector{String}
    msp_files = readdlm(input4) |> vec |> Vector{String}
else
    error(
        "Problem with input. Input 1 and 2 should both be VCF files (ends with .vcf or .vcf.gz), " * 
        "Or they should both point to a file which contains a list of VCF files."
    )
end
df, _ = get_ancestry_specific_r2(
    genotype_files, imputed_files, summary_files, msp_files, input5,
    ancestry_names=ancestry_names
)

# save as dataframe
outdir = dirname(outfile)
isdir(outdir) || mkpath(outdir)
CSV.write(outfile, df, delim='\t')
println("Done! Ancestry-specific aggregate R2 saved to $outfile")
