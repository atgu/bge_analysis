# usage: julia concordance.jl --truth file1 --impt file2 --out file3
# --truth: Full file path to the ground truth genotypes (must end in .vcf or .vcf.gz), 
#         or a file containing a list of VCF files
# --impt: Full file path to the imputed genotypes (must end in .vcf or .vcf.gz), 
#         or a file containing a list of VCF files
# --out: Full file path to output file name. Output will be non-reference
#         concordance (i.e the percent of non-reference alleles that were
#         correctly imputed)
# --summary: An optional input. Either Full file path to a comma or tab separated
#         summary file (header must include at least 6 columns with the names CHR/POS/REF/ALT/AF/isImputed), 
#         or a file containing a list of summary files
# --maf-bins: Comma-separated list of minor allele frequencies used to bin SNPs.
#             Defaults to 0.0,0.0005,0.001,0.004,0.0075,0.0125,0.04,0.1,0.2,0.5

include("utilities.jl")

# parse inputs to this script
using ArgParse
s = ArgParseSettings()
@add_arg_table! s begin
    "--truth"
        required = true
        arg_type = String
    "--impt"
        required = true
        arg_type = String
    "--out"
        required = true
        arg_type = String
    "--summary"
        default = ""
        arg_type = String
    "--maf-bins"
        default = "0.0,0.0005,0.001,0.004,0.0075,0.0125,0.04,0.1,0.2,0.5"
        arg_type = String
end
parsed_args = parse_args(s)
input1 = parsed_args["truth"]
input2 = parsed_args["impt"]
outfile = parsed_args["out"]
maf_bins = split(parsed_args["maf-bins"], ',')

# check if input is a list of files
input1_isvcf = endswith(input1, ".vcf") || endswith(input1, ".vcf.gz")
input2_isvcf = endswith(input2, ".vcf") || endswith(input2, ".vcf.gz")
if input1_isvcf && input2_isvcf
    genotype_file = input1
    imputed_file = input2
    summary_file = parsed_args["summary"]
elseif !input1_isvcf && !input2_isvcf # input is list of files
    genotype_file = readdlm(input1) |> vec |> Vector{String}
    imputed_file = readdlm(input2) |> vec |> Vector{String}
    summary_file = parsed_args["summary"] == "" ? 
        ["" for _ in eachindex(genotype_file)] : 
        (readdlm(summary_file) |> vec |> Vector{String})
else
    error(
        "Problem with input. Input 1 and 2 should both be VCF files (ends " * 
        "with .vcf or .vcf.gz), or they should both point to a file which " * 
        "contains a list of VCF files."
    )
end

@info "Computing non-reference concordance"
concordances, counters, mafs = compute_non_ref_concordance(
    genotype_file, imputed_file, summary_file=summary_file
)

# save as dataframe
df = DataFrame("concordance" => concordances, 
               "Num_true_positives" => counters,
               "mafs" => mafs
)
outdir = dirname(outfile)
isdir(outdir) || mkpath(outdir)
CSV.write(outfile, df, delim='\t')
println("Done! Concordance file saved to $outfile")
