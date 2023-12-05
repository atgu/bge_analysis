using VCFTools
using DataFrames
using CSV
using DelimitedFiles
using Statistics
using ProgressMeter

"""
The SNP struct defines a "SNP" supporting equality comparisons. We can compare
if 2 SNP is the same by checking chr, pos, ref, and alt. For 2 SNPs to be 
considered equal, they must have the same chr and pos, but if ref/alt is
flipped, we still consider them the "same SNP". If we have a `Vector{SNP}`, 
then `indexin` and `intersect` are overloaded from Julia Base. 

Note: if a record is multiallelic, then there will be multiple alt alleles. One
can still use the SNP struct by combining the "sorted" alt alleles via ','. For 
example, if a SNP have ALT alleles 'T' and 'A', then the `alt == "A,T"`
"""
struct SNP
    chr::String
    pos::Int
    ref::String
    alt::String
end
# more complex constructors
function SNP(chr::Int, pos::Int, ref::AbstractString, alt::AbstractVector)
    return SNP(string(chr), pos, ref, join(sort!(alt),','))
end
function SNP(chr::AbstractString, pos::Int, ref::AbstractString, alt::AbstractVector)
    return SNP(chr, pos, ref, join(sort!(alt),','))
end
function SNPs(chrs::AbstractVector, poss::AbstractVector{Int}, 
    refs::AbstractVector, alts::AbstractVector)
    snps = SNP[]
    for (chr, pos, ref, alt) in zip(chrs, poss, refs, alts)
        push!(snps, SNP(chr, pos, ref, alt))
    end
    return snps
end

# methods supported by SNP struct
import Base.==, Base.indexin, Base.intersect
function ==(x::SNP, y::SNP)
    alleles_match = ((x.ref == y.ref) && (x.alt == y.alt)) || 
        ((x.ref == y.alt) && (x.alt == y.ref))
    if (x.chr != y.chr) || (x.pos != y.pos) || !alleles_match
        return false
    end
    return true
end
function indexin(a::AbstractArray{SNP}, b::AbstractArray{SNP})
    # same as Base.index but works on a vector of SNP.
    # Note to self: it is important to include both SNP(chr, pos, ref, alt) and SNP(chr, pos, alt, ref)
    # when checking if a SNP exist in the other vector
    inds = keys(b)
    bdict = Dict{eltype(b),eltype(inds)}()
    for (snp, ind) in zip(b, inds)
        get!(bdict, snp, ind)
        get!(bdict, SNP(snp.chr, snp.pos, snp.alt, snp.ref), ind)
    end
    result = Union{eltype(inds), Nothing}[]
    for snp in a
        result1 = get(bdict, snp, nothing)
        result2 = get(bdict, SNP(snp.chr, snp.pos, snp.alt, snp.ref), nothing)
        if !isnothing(result1)
            push!(result, result1)
        elseif !isnothing(result2)
            push!(result, result2)
        else
            push!(result, nothing)
        end
    end
    return result
end
function intersect(a::AbstractArray{SNP}, b::AbstractArray{SNP})
    c = SNP[]
    for snp in a
        snp in b && push!(c, snp)
    end
    return unique!(c)
end

"""
    get_AF(x::AbstractMatrix)

Compute alternate allele frequency on dosage matrix `x`, 
skipping missing values if there are any
"""
function get_AF(x::AbstractMatrix)
    return mean.(skipmissing.(eachcol(x))) ./ 2
end

"""
    get_ancestry_names(mspfile::AbstractString)

Reads the super population names from output of `Rfmix2`. Assumes the first line
takes the form `#Subpopulation order/codes: POP1=0\tPOP2=1\tPOP3=2` for 3 populations
"""
function get_ancestry_names(mspfile::AbstractString)
    l = readline(mspfile)
    codes = split(l[29:end], '\t')
    ancestries = [split(code, '=')[1] for code in codes] |> Vector{String}
    return ancestries
end

"""
    import_Xtrue_Ximp(genotype_file, imputed_file, summary_file; [use_dosage])

Imports VCF files `genotype_file` and `imputed_file` into numeric matrices and
subset them so their are on the same set of SNPs. If REF/ALT is opposite in 
`genotype_file` and `imputed_file`, we will flip corresponding entries in 
`genotype_file` so that the dosages from `genotype_file` and genotypes of 
`imputed_file` can be correlated directly. Minor allele frequency is computed
based on `imputed_file`
"""
function import_Xtrue_Ximp(
    genotype_file::String, # VCF (will import GT field)
    imputed_file::String; # VCF (import DS = dosage field, unless use_dosage=false, in which case we import GT)
    summary_file::String = "", # must contain CHR/POS/REF/ALT/isImputed columns
    use_dosage::Bool=true
    )
    # import array genotype data
    array_data, array_sampleID, array_chr, array_pos, _, array_ref, array_alt = 
        convert_gt(Float64, genotype_file, save_snp_info=true, 
        msg="Importing array genotype VCF data")
    array_snps = SNPs(array_chr, array_pos, array_ref, array_alt)

    # import imputed data
    import_func = use_dosage ? convert_ds : convert_gt
    imputed_data, imputed_sampleID, imputed_chr, imputed_pos, _, imputed_ref, imputed_alt = 
        import_func(Float64, imputed_file, save_snp_info=true, 
        msg="Importing imputed VCF")
    imputed_snps = SNPs(imputed_chr, imputed_pos, imputed_ref, imputed_alt)
    
    # compute maf (between 0 and 0.5) from imputed data
    AF = get_AF(imputed_data)
    MAF = copy(AF)
    idx = findall(x -> x > 0.5, MAF)
    MAF[idx] .-= 1 .- MAF[idx]

    # shared samples
    shared_samples = intersect(imputed_sampleID, array_sampleID)
    println("matched $(length(shared_samples)) samples")

    # shared SNPs
    shared_snps = intersect(imputed_snps, array_snps)
    if summary_file != ""
        summary_df = CSV.read(summary_file, DataFrame)
        keep_idx = findall(x -> x == false, summary_df[!, "isImputed"])
        keep_chr = summary_df[keep_idx, "CHR"]
        keep_pos = summary_df[keep_idx, "POS"]
        keep_ref = summary_df[keep_idx, "REF"]
        keep_alt = summary_df[keep_idx, "ALT"]
        keep_snps = SNPs(keep_chr, keep_pos, keep_ref, keep_alt)
        idx = findall(x -> !(x in keep_snps), shared_snps)
        shared_snps = shared_snps[idx]
    end
    println("matched $(length(shared_snps)) SNPs")

    # subset array/imputed data so only matching SNPs remain
    array_col_idx = indexin(shared_snps, array_snps)
    imputed_col_idx = indexin(shared_snps, imputed_snps)
    array_row_idx = indexin(shared_samples, array_sampleID)
    imputed_row_idx = indexin(shared_samples, imputed_sampleID)
    Xgeno = array_data[array_row_idx, array_col_idx]
    Ximpt = imputed_data[imputed_row_idx, imputed_col_idx]
    size(Ximpt, 2) == size(Xgeno, 2) == length(shared_snps) || 
        error("expected size(Ximpt, 2) == size(Xgeno, 2) == length(shared_snps)")

    # flip 2 to 0 and 0 to 2 if ref/alt allele is opposite
    imputed_snps = imputed_snps[imputed_col_idx]
    imputed_ref = imputed_ref[imputed_col_idx]
    imputed_alt = imputed_alt[imputed_col_idx]
    MAF = MAF[imputed_col_idx]
    array_snps = array_snps[array_col_idx]
    for (i, (imputed_snp, genotyped_snp)) in enumerate(zip(imputed_snps, array_snps))
        if (imputed_snp.ref == genotyped_snp.ref) && (imputed_snp.alt == genotyped_snp.alt)
            continue
        elseif (imputed_snp.ref == genotyped_snp.alt) && (imputed_snp.alt == genotyped_snp.ref) # flip
            Xgeno[:, i] .= 2.0 .- Xgeno[:, i]
        else
            error("shouldn't happen!")
        end
    end

    size(Xgeno) == size(Ximpt) || error("expected size(Xgeno) == size(Ximpt)")
    length(MAF) == length(shared_snps) || error("expected length(MAF) == length(shared_snps)")
    size(Xgeno, 2) == length(shared_snps) || error("expected size(Xgeno, 2) == length(shared_snps)")

    return Xgeno, Ximpt, MAF, shared_samples, shared_snps
end

function get_aggregate_R2(
    genotype_file::String,
    imputed_file::String;
    summary_file::String = "", # must contain CHR/POS/REF/ALT/isImputed columns
    maf_bins::Vector{Float64}=[0.0, 0.0005, 0.001, 0.004, 0.0075, 0.0125, 0.04, 0.1, 0.2, 0.5]
    )
    # import data
    Xgeno, Ximpt, mafs, _, _ = import_Xtrue_Ximp(
        genotype_file, imputed_file; summary_file=summary_file
    )

    # compute aggregate R2
    R2 = Float64[]
    nsnps = Int[]
    for i in 2:length(maf_bins)
        # find SNPs in current bin
        bin_low, bin_high = maf_bins[i-1], maf_bins[i]
        idx = findall(x -> bin_low ≤ x < bin_high, mafs)
        if length(idx) == 0
            @warn("maf_bin [$bin_low, $bin_high] has no SNP! Assigning R2 of NaN")
            push!(R2, NaN)
            push!(nsnps, 0)
            continue
        end

        # compute squared pearson correlation
        truth = vec(Xgeno[:, idx])
        imptd = vec(Ximpt[:, idx])
        non_missing_idx = intersect(findall(!ismissing, truth), findall(!ismissing, imptd))
        my_R2 = cor(truth[non_missing_idx], imptd[non_missing_idx])^2
        push!(R2, my_R2)
        push!(nsnps, length(idx))
    end
    return R2, nsnps
end

function get_aggregate_R2(
    genotype_file::Vector{String},
    imputed_file::Vector{String};
    summary_file::Vector{String} = ["" for _ in eachindex(genotype_files)], # must contain CHR/POS/REF/ALT/isImputed columns
    maf_bins::Vector{Float64}=[0.0, 0.0005, 0.001, 0.004, 0.0075, 0.0125, 0.04, 0.1, 0.2, 0.5]
    )
    length(genotype_file) == length(imputed_file) == length(summary_file) || 
        error(
            "Expected genotype_file, imputed_file, summary_file to have the same length"
        )

    # compute aggregate R2 over all files
    agg_R2 = zeros(length(maf_bins) - 1)
    tot_snps = zeros(Int, length(maf_bins) - 1)
    @showprogress for (g, i, s) in zip(genotype_file, imputed_file, summary_file)
        R2, nsnps = get_aggregate_R2(
            g, i, summary_file=s, maf_bins=maf_bins
        )
        agg_R2 .+= R2 .* nsnps
        tot_snps .+= nsnps
    end
    agg_R2 ./= tot_snps

    return agg_R2, tot_snps
end


"""
    create_LAI_mapping_matrices(mspfile::String, n_ancestries::Int)

Helper function for creating masking matrices. 

# Inputs 
+ `mspfile`: msp file outputted from running rfmix2, see Tractor tutorial:
    https://github.com/Atkinson-Lab/Tractor-tutorial/blob/main/Rfmix.md
+ `n_ancestries`: Number of ancestries used in rfmix2
"""
function create_LAI_mapping_matrices(
    mspfile::String, 
    n_ancestries::Int;
    ancestry_names::Vector{String} = ["$(ancestry)$(ancestry)" for ancestry in 0:n_ancestries-1]
    )
    length(ancestry_names) == n_ancestries || 
        error("expected length(ancestry_names) == n_ancestries")
    msp = CSV.read(mspfile, DataFrame, header=2)
    mat = msp[:, 7:end]
    segment_nsnps = msp[!, "n snps"]
    nsnps = sum(segment_nsnps)
    nsamples = size(mat, 2) >> 1
    per_snp_ancestry_background = Dict{String, BitMatrix}()

    for ancestry1 in 0:n_ancestries-1, ancestry2 in ancestry1:n_ancestries-1
        bm = falses(nsamples, nsnps)
        # (i, j) indexes into the n by p bitmatrix `bm`.
        # `bm[i, j] = true` if the genotype for sample i at SNP j matches
        # ancestry1 and ancestry2, otherwise `bm[i, j] = 0`. E.g. if
        # (ancestry1, ancestry2) = (0, 0) = (AFR, AFR) and (i, j) = (1, 1)
        # then `bm[i, j] = 1` if have both haplotypes for sample `i` at SNP `j`
        # have AFR background
        for i in 1:nsamples, jj in eachindex(segment_nsnps)
            # jj indexes into the rows of the msp matrix
            # recall each row of msp matrix is a "ancestry segment" containing 
            # a number of SNPs all of which have the same ancestry background
            anc1 = mat[jj, 2i - 1]
            anc2 = mat[jj, 2i]
            if (ancestry1 == anc1 && ancestry2 == anc2) || (ancestry1 == anc2 && ancestry2 == anc1)
                j_start = jj == 1 ? 1 : sum(segment_nsnps[1:jj-1]) + 1
                j_end = sum(segment_nsnps[1:jj])
                for j in j_start:j_end
                    bm[i, j] = true
                end
            end
        end
        name = ancestry_names[ancestry1+1] * '_' * ancestry_names[ancestry2+1]
        per_snp_ancestry_background[name] = bm
    end

    # read the sample IDs for masking matrix
    sampleID = String[]
    for name in names(msp)[7:end]
        push!(sampleID, name[1:end-2]) # remove last ".0" or ".1" 
    end
    unique!(sampleID)    

    return per_snp_ancestry_background, sampleID, nsnps
end

function get_ancestry_specific_r2(
    genotype_file::String, # VCF (will import GT field)
    imputed_file::String, # VCF (import DS field, unless use_dosage=false, in which case we import GT)
    mspfile::String, # output of Rfmix2 (input to Rfmix2 MUST be genotype_file)
    ancestry_names::Vector{String};
    summary_file::String = "", # must contain CHR/POS/REF/ALT/AF/isImputed columns
    use_dosage::Bool=true, # whether to import imputed data as dosage or genotypes
    maf_bins::Vector{Float64}=[0.0, 0.0005, 0.001, 0.004, 0.0075, 0.0125, 0.04, 0.1, 0.2, 0.5],
    )
    # compute background ancestries on the phased genotypes
    n_ancestries = length(ancestry_names)
    ancestry_masks, _, nsnps = create_LAI_mapping_matrices(
        mspfile, n_ancestries, ancestry_names=ancestry_names)

    # before continueing, check for a weird error I got on chr10 of PAISA data
    if nsnps != nrecords(genotype_file)
        error(
            "The `n snps` column in msp file does not sum to the number " * 
            "of SNPs in $genotype_file. I got this weird error before, so " * 
            "this might be Rfmix2 filtering out SNPs based on unknown criteria"
        )
    end

    # import the genotype and imputed data matrices, after matching them
    Xgeno, Ximpt, mafs, samples, snps = import_Xtrue_Ximp(
        genotype_file, imputed_file, summary_file=summary_file, 
        use_dosage=use_dosage, 
    )

    # Note: background ancestries were computed on the phased genotypes, 
    # which gives a 0/1 bit indicating whether the given sample/SNP has a 
    # specific background ancestry. However, the masking matrices in 
    # `ancestry_masks` does not come with sample/SNP labeling, so we must 
    # import that information separately
    _, array_sampleIDs, array_chr, array_pos, _, array_ref, array_alt = 
        convert_gt(Float64, genotype_file, save_snp_info=true, 
        msg="Importing VCF data")
    array_snps = SNPs(array_chr, array_pos, array_ref, array_alt)

    # figure out which sample/SNPs can be used and subset the masks
    row_idx = indexin(samples, array_sampleIDs)
    col_idx = indexin(snps, array_snps)
    for (ancestry, has_this_ancestry) in ancestry_masks
        ancestry_masks[ancestry] = has_this_ancestry[row_idx, col_idx]
    end

    # prepare dataframe
    df = DataFrame(maf_bins = String[])
    for i in 1:length(maf_bins)-1
        if i == length(maf_bins)
            push!(df[!, "maf_bins"], "[$(maf_bins[i]), $(maf_bins[i+1])]")
        else
            push!(df[!, "maf_bins"], "[$(maf_bins[i]), $(maf_bins[i+1]))")
        end
    end

    # loop through each background ancestry
    nsnps = zeros(length(maf_bins) - 1, length(ancestry_masks))
    counter = 1
    for (ancestry, has_this_ancestry) in ancestry_masks
        # compute aggregate R2
        R2 = Float64[]
        for i in 2:length(maf_bins)
            # find SNPs in current bin
            bin_low, bin_high = maf_bins[i-1], maf_bins[i]
            idx = findall(x -> bin_low ≤ x < bin_high, mafs)
            if length(idx) == 0
                @warn("maf_bin [$bin_low, $bin_high] has no SNP! Assining R2 of NaN")
                push!(R2, NaN)
                continue
            end

            # compute squared pearson correlation
            mask = has_this_ancestry[:, idx]
            truth = vec(Xgeno[:, idx][mask])
            imptd = vec(Ximpt[:, idx][mask])
            non_missing_idx = intersect(findall(!ismissing, truth), findall(!ismissing, imptd))
            push!(R2, cor(truth[non_missing_idx], imptd[non_missing_idx])^2)
            nsnps[i-1, counter] = length(idx) # record number of SNPs used to compute R2
        end
        df[!, ancestry] = R2
        counter += 1
    end

    # returns a dataframe containing aggregate R2 for different ancestries
    # the second return argument is the number of SNPs in each bin used to 
    # calculate the given R2
    return df, nsnps
end

# computes aggregate R2, averaging over all files (but taking into account the 
# number of SNPs in each MAF bin). Thus, this is not a "true aggregate R2" over
# all files. 
function get_ancestry_specific_r2(
    genotype_file::Vector{String}, # VCF files (will import GT field), e.g. one for each chr
    imputed_file::Vector{String}, # VCF files (import DS field, unless use_dosage=false, in which case we import GT)
    msp_file::Vector{String}, # output of Rfmix2 (input to Rfmix2 MUST be genotype_file)
    ancestry_names::Vector{String};
    summary_file::Vector{String} = ["" for _ in eachindex(genotype_file)], 
    use_dosage::Bool=true,
    maf_bins::Vector{Float64}=[0.0, 0.0005, 0.001, 0.004, 0.0075, 0.0125, 0.04, 0.1, 0.2, 0.5],
    )
    # get ancestry specific R2 for first set of files
    df, nsnps = get_ancestry_specific_r2(
        genotype_file[1], imputed_file[1], msp_file[1], ancestry_names, 
        summary_file = summary_file[1], use_dosage=use_dosage, 
        maf_bins=maf_bins
    )
    R2s = df[:, 2:end] |> Matrix{Float64}

    # if a certain ancestry/maf-bin have 0 SNPs, R2 will be NaN. Set those to 0
    idx = findall(isnan, R2s)
    R2s[idx] .= 0
    nsnps[idx] .= 0

    # variable to return: average R2 over all files, accounting bin sizes
    tot_snps = copy(nsnps)
    avg_R2 = copy(R2s) .* tot_snps

    # loop over remaining file
    for i in 2:length(genotype_file)
        df, nsnps = get_ancestry_specific_r2(
            genotype_file[i], imputed_file[i], msp_file[i], ancestry_names, 
            summary_file = summary_file[i], use_dosage=use_dosage, 
            maf_bins=maf_bins, 
        )
        R2s = df[!, 2:end] |> Matrix{Float64}
        idx = findall(!isnan, R2s)

        avg_R2[idx] .+= R2s[idx] .* nsnps[idx]
        tot_snps[idx] .+= nsnps[idx]
    end

    # final dataframe
    avg_R2 ./= tot_snps
    df = copy(df)
    df[:, 2:end] .= avg_R2

    return df, tot_snps
end

function compute_precision_sensitivity(ximp::AbstractMatrix, xtrue::AbstractMatrix)
    size(ximp) == size(xtrue) || error("check dimension")
    n, p = size(ximp)
    precisions, sensitivities = Float64[], Float64[]
    for j in 1:p
        TP, FP, FN = 0, 0, 0
        for i in 1:n
            ismissing(xtrue[i, j]) && continue
            ismissing(ximp[i, j]) && continue

            # true positives
            if xtrue[i, j] == ximp[i, j] == 1
                TP += 1
            elseif xtrue[i, j] == 1 && ximp[i, j] == 2
                TP += 1
            elseif xtrue[i, j] == 2 && ximp[i, j] == 1
                TP += 1
            elseif xtrue[i, j] == ximp[i, j] == 2
                TP += 2
            end

            # false positives
            if xtrue[i, j] == 0 && ximp[i, j] == 1
                FP += 1
            elseif xtrue[i, j] == 0 && ximp[i, j] == 2
                FP += 2
            elseif xtrue[i, j] == 1 && ximp[i, j] == 2
                FP += 1
            end

            # false negatives
            if xtrue[i, j] == 1 && ximp[i, j] == 0
                FN += 1
            elseif xtrue[i, j] == 2 && ximp[i, j] == 1
                FN += 1
            elseif xtrue[i, j] == 2 && ximp[i, j] == 0
                FN += 2
            end
        end
        push!(precisions, TP / (TP + FP))
        push!(sensitivities, TP / (TP + FN))
    end
    return precisions, sensitivities
end

function compute_non_ref_concordance(ximp::AbstractMatrix, xtrue::AbstractMatrix)
    size(ximp) == size(xtrue) || error("check dimension")
    n, p = size(ximp)
    concordance = zeros(p)
    counters = zeros(Int, p)
    for j in 1:p
        tot_nonref_alleles = sum(skipmissing(@view(xtrue[:, j])))
        TP = 0
        for i in 1:n
            ismissing(xtrue[i, j]) && continue
            ismissing(ximp[i, j]) && continue

            # true positives
            if xtrue[i, j] == ximp[i, j] == 1
                TP += 1
            elseif xtrue[i, j] == 1 && ximp[i, j] == 2
                TP += 1
            elseif xtrue[i, j] == 2 && ximp[i, j] == 1
                TP += 1
            elseif xtrue[i, j] == ximp[i, j] == 2
                TP += 2
            end
        end
        concordance[j] = TP / tot_nonref_alleles
        counters[j] = TP
    end
    return concordance, counters
end

function compute_non_ref_concordance(
    genotype_file::String, # VCF (will import GT field)
    imputed_file::String; # VCF (import DS field, unless use_dosage=false, in which case we import GT)
    summary_file::String = ""
    )

    # import the full genotype and imputed data matrices
    Xgeno, Ximpt, mafs, _, _ = import_Xtrue_Ximp(
        genotype_file, imputed_file, summary_file=summary_file, use_dosage=false,
    )

    concordance, counters = compute_non_ref_concordance(Ximpt, Xgeno)

    return concordance, counters, mafs
end

function compute_non_ref_concordance(
    genotype_file::Vector{String}, # VCF (will import GT field)
    imputed_file::Vector{String}; # VCF (import DS field, unless use_dosage=false, in which case we import GT)
    summary_file::Vector{String} = ["" for _ in eachindex(genotype_file)]
    )
    concordances, counters, mafs = Float64[], Int[], Float64[]
    for (gtfile, impfile, summ_file) in zip(genotype_file, imputed_file, summary_file)
        # import the full genotype and imputed data matrices
        Xgeno, Ximpt, mf, _, _ = import_Xtrue_Ximp(
            gtfile, impfile, summary_file=summ_file, use_dosage=false,
        )

        # compute confordance for current file and save 
        conc, cont = compute_non_ref_concordance(Ximpt, Xgeno)
        append!(concordances, conc)
        append!(counters, cont)
        append!(mafs, mf)
        length(concordances) == length(counters) == length(mafs) || 
            error("length different!")
    end

    return concordances, counters, mafs
end
