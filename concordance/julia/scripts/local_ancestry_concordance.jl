# usage: julia local_ancestry_concordance.jl --truth file1 --impt file2 --msp file3 --out file4
# --truth: Full file path to the ground truth genotypes (must end in .vcf or .vcf.gz)
# --impt: Full file path to the imputed genotypes (must end in .vcf or .vcf.gz)
# --msp: Full file path to the .msp.tsv file (output of rfmix2) 
# --out: Full file path to output file name (no extensions). Each ancestry will 
#        output a separate dataframe including the concordances for each SNP. 
# --summary: An optional which allows users to NOT consider certain SNPs that 
#            exist in ground truth or imputed genotype files. Must be a
#            full file path to a comma or tab separated summary file (header must
#            include at least 5 columns with the names CHR/POS/REF/ALT/isImputed)
#            The scripts are set up so that only SNPs listed as `true` in the 
#            `isImputed` column will be considered. 

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
    "--msp"
        required = true
        arg_type = String
    "--out"
        required = true
        arg_type = String
    "--summary"
        default = ""
        arg_type = String
end
parsed_args = parse_args(s)
input1 = parsed_args["truth"]
input2 = parsed_args["impt"]
input3 = parsed_args["msp"]
outfile = parsed_args["out"]

# check if input is a list of files
input1_isvcf = endswith(input1, ".vcf") || endswith(input1, ".vcf.gz")
input2_isvcf = endswith(input2, ".vcf") || endswith(input2, ".vcf.gz")
if input1_isvcf && input2_isvcf
    genotype_file = input1
    imputed_file = input2
    msp_file = input3
    summary_file = parsed_args["summary"]
    ancestry_names = get_ancestry_names(msp_file)
else
    error(
        "Problem with input. Input 1 and 2 should both be VCF files (ends " * 
        "with .vcf or .vcf.gz). Note this script does not compute over a " * 
        "a list of VCF files. One should submit a separate job for each chromosome."
    )
end

@info "Computing concordance by local ancestries"
result = get_ancestry_specific_concordance(
    genotype_file, imputed_file, msp_file, ancestry_names=ancestry_names,
    summary_file=summary_file
)

# save each dataframe (one for each ancestry) as dataframe
outdir = dirname(outfile)
isdir(outdir) || mkpath(outdir)
for (ancestry, df) in result
    CSV.write(outfile * ancestry, df, delim='\t')
end
println("Done! Ancestry-specific concordances saved to $outdir")
