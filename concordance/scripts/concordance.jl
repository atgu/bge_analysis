# usage: julia concordance.jl input1 input2 input3 input4
# input1: Full file path to the ground truth genotypes (must end in .vcf or .vcf.gz), 
#         or a file containing a list of VCF files
# input2: Full file path to the imputed genotypes (must end in .vcf or .vcf.gz), 
#         or a file containing a list of VCF files
# input3: Full file path to a comma or tab separated summary file (header must 
#         include at least 6 columns with the names CHR/POS/REF/ALT/AF/isImputed), 
#         or a file containing a list of summary files
# input4: Full file path to output file name. Output will be non-reference
#         concordance (i.e the percent of non-reference alleles that were
#         correctly imputed)

include("utilities.jl")

# inputs to this script
input1 = ARGS[1]
input2 = ARGS[2]
input3 = ARGS[3]
outfile = ARGS[4]

# default minor allele frequency bins
maf_bins = [0.0, 0.0005, 0.00100, 0.00400, 0.0075, 0.0125, 0.04, 0.1, 0.2, 0.5]

# compute non-reference concordance
@info "Computing non-reference concordance"
input1_isvcf = endswith(input1, ".vcf") || endswith(input1, ".vcf.gz")
input2_isvcf = endswith(input2, ".vcf") || endswith(input2, ".vcf.gz")
if input1_isvcf && input2_isvcf
    genotype_files = input1
    imputed_files = input2
    summary_files = input3
elseif !input1_isvcf && !input2_isvcf
    genotype_files = readdlm(input1) |> vec |> Vector{String}
    imputed_files = readdlm(input2) |> vec |> Vector{String}
    summary_files = readdlm(input3) |> vec |> Vector{String}
else
    error(
        "Problem with input. Input 1 and 2 should both be VCF files (ends with .vcf or .vcf.gz), " * 
        "Or they should both point to a file which contains a list of VCF files."
    )
end
concordances, counters, mafs = compute_non_ref_concordance(
    genotype_files, imputed_files, summary_files
)

# save as dataframe
df = DataFrame("concordance" => concordances, 
               "number of true positives" => counters,
               "mafs" => mafs
)
outdir = dirname(outfile)
isdir(outdir) || mkpath(outdir)
CSV.write(outfile, df, delim='\t')
println("Done! Concordance file saved to $outfile")
