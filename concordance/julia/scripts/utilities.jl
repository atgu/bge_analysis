using VCFTools
using DataFrames
using CSV
using DelimitedFiles
using Statistics
using ProgressMeter

include("struct.jl")

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
    import_Xtrue_Ximp(genotype_file, imputed_file; [summary_file], [use_dosage])

Imports VCF files `genotype_file` and `imputed_file` into numeric matrices and
subset them so their are on the same set of SNPs. If REF/ALT is opposite in 
`genotype_file` and `imputed_file`, we will flip corresponding entries in 
`genotype_file` so that the dosages from `genotype_file` and genotypes of 
`imputed_file` can be correlated directly. Minor allele frequency is computed
    based on `genotype_file`
"""
function import_Xtrue_Ximp(
    genotype_file::String, # VCF (will import GT field)
    imputed_file::String; # VCF (import DS = dosage field, unless use_dosage=false, in which case we import GT)
    summary_file::String = "", # must contain CHR/POS/REF/ALT/isImputed columns
    use_dosage::Bool=true,
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
    AF = get_AF(array_data)
    MAF = copy(AF)
    idx = findall(x -> x > 0.5, MAF)
    MAF[idx] .= 1 .- MAF[idx]

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
    array_snps = array_snps[array_col_idx]
    MAF = MAF[array_col_idx]
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

function aggregate_R2(Xgeno::AbstractVecOrMat, Ximpt::AbstractVecOrMat)
    truth = vec(Xgeno)
    imptd = vec(Ximpt)
    idx = intersect(findall(!ismissing, truth), findall(!ismissing, imptd))
    return cor(truth[idx], imptd[idx])^2
end

function get_aggregate_R2(
    genotype_file::String,
    imputed_file::String;
    summary_file::String = "", # must contain CHR/POS/REF/ALT/isImputed columns
    maf_bins::Vector{Float64}=[0.0, 0.0005, 0.001, 0.005, 0.01, 0.02, 0.03, 0.04, 0.05, 0.1, 0.2, 0.3, 0.4, 0.5]
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
        my_R2 = aggregate_R2(Xgeno[:, idx], Ximpt[:, idx])
        push!(R2, my_R2)
        push!(nsnps, length(idx))
    end
    return R2, nsnps
end

function get_aggregate_R2(
    genotype_file::Vector{String},
    imputed_file::Vector{String};
    summary_file::Vector{String} = ["" for _ in eachindex(genotype_files)], # must contain CHR/POS/REF/ALT/isImputed columns
    maf_bins::Vector{Float64}=[0.0, 0.0005, 0.001, 0.005, 0.01, 0.02, 0.03, 0.04, 0.05, 0.1, 0.2, 0.3, 0.4, 0.5]
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
    genotyped_samples::Vector{String},
    genotyped_snps::Vector{SNP},
    mspfile::String, 
    n_ancestries::Int;
    ancestry_names::Vector{String} = ["$(ancestry)$(ancestry)" for ancestry in 0:n_ancestries-1]
    )
    # include("/u/home/b/biona001/bge_analysis/concordance/scripts/utilities.jl")
    # genotype_file = "/u/home/b/biona001/project-loes/gpc_gsa/array/all_GPC_cohorts.chr20.sharedSNPs.vcf.gz"
    # imputed_file = "/u/home/b/biona001/project-loes/gpc_gsa/imputed/typedSNPs/chr20.vcf.gz"
    # mspfile = "/u/home/b/biona001/project-loes/gpc_gsa/tbb/gpc_gsa_w2.chr20.typed.deconvolved.msp.tsv"
    # genotype_file = "/u/home/b/biona001/project-loes/ForBen_genotypes_subset/QC_hg38_conformed_king/chr22.vcf.gz"
    # imputed_file = "/u/home/b/biona001/project-loes/GLIMPSE2_toni/typedSNPs/chr22.sampleQC.snpQC.vcf.gz"
    # mspfile = "/u/home/b/biona001/project-loes/ForBen_genotypes_subset/LAI/output_v2/chr22.msp.tsv"
    # summary_file = ""
    # ancestry_names = get_ancestry_names(mspfile)
    # use_dosage = true
    # maf_bins = collect(LinRange(0, 0.5, 100))
    # n_ancestries = length(ancestry_names)
    # Xgeno, Ximpt, mafs, genotyped_samples, genotyped_snps = import_Xtrue_Ximp(
    #     genotype_file, imputed_file, summary_file=summary_file, 
    #     use_dosage=use_dosage, 
    # )
    # @time per_snp_ancestry_background = create_LAI_mapping_matrices(
    #     genotyped_samples, genotyped_snps, mspfile, n_ancestries, 
    #     ancestry_names=ancestry_names
    #     ) # 189 seconds (new code = 7.58 seconds)
    # @test findall(per_snp_ancestry_background["EUR_EUR"][1:20]) == [1, 2, 3, 4, 6, 7, 10, 15, 17, 19]
    # @test findall(per_snp_ancestry_background["EUR_EUR"][2341:2380]) == [1, 3, 7, 8, 9, 12, 13, 14, 20, 23, 25, 26, 28, 30, 31, 33, 34, 40]
    # @test findall(per_snp_ancestry_background["EUR_EUR"][end-20:end]) == [1, 2, 3, 6, 9, 10, 11, 14, 16, 18, 21]
    # @test findall(per_snp_ancestry_background["AFR_EUR"][1:300]) == [8, 9, 12, 14, 18, 33, 59, 79, 91, 117, 118, 125, 156, 162, 200, 242, 252, 255, 276, 281, 288, 289]

    length(ancestry_names) == n_ancestries || 
        error("expected length(ancestry_names) == n_ancestries")
    msp = CSV.read(mspfile, DataFrame, header=2)
    mat = msp[:, 7:end]
    spos = msp[!, "spos"] # start pos of each LAI segment
    epos = msp[!, "epos"] # end pos of each LAI segment
    nsnps = length(genotyped_snps)
    nsamples = length(genotyped_samples)
    msp_samples = split.(names(msp)[7:end], '.')
    msp_samples = [x[1] for x in msp_samples] |> Vector{String}
    genotyped_pos = [snp.pos for snp in genotyped_snps]
    per_snp_ancestry_background = Dict{String, BitMatrix}()
    bm2mat_sampleidx = indexin(genotyped_samples, msp_samples)

    # check for errors
    all(!isnothing, bm2mat_sampleidx) || 
        error("some sample ID in array data does not have LAI computed in msp file")
    issorted(spos) && issorted(epos) || error("spos in mspfile $mspfile is not sorted!")

    for ancestry1 in 0:n_ancestries-1, ancestry2 in ancestry1:n_ancestries-1
        bm = falses(nsamples, nsnps)
        # (i, j) indexes into the n by p bitmatrix `bm`
        # `bm[i, j] = true` if the genotype for sample `genotyped_samples[i]` at
        # SNP `genotyped_snps[j]` matches ancestry1 and ancestry2, 
        # otherwise `bm[i, j] = 0`. For example, if
        # (ancestry1, ancestry2) = (0, 0) = (AFR, AFR) and (i, j) = (1, 1)
        # then `bm[i, j] = 1` if have both haplotypes for sample 
        # `genotyped_samples[i]` at SNP `genotyped_snps[j]` have AFR background
        @inbounds @showprogress for jj in eachindex(spos)
            # (jj, ii) indexes into `mat` (the msp matrix)
            # recall each row of msp matrix is a "ancestry segment" containing 
            # a number of SNPs all of which have the same ancestry background
            segment = spos[jj]:(epos[jj] - 1)
            for j in findall(x -> x ∈ segment, genotyped_pos)
                for (i, ii) in zip(1:nsamples, bm2mat_sampleidx)
                    anc1 = mat[jj, ii]
                    anc2 = mat[jj, ii + 1]
                    anc_match = (ancestry1 == anc1) && (ancestry2 == anc2)
                    anc_cross_match = (ancestry1 == anc2) && (ancestry2 == anc1)
                    if anc_match || anc_cross_match
                        bm[i, j] = true
                    end
                end
            end
        end

        name = ancestry_names[ancestry1+1] * '_' * ancestry_names[ancestry2+1]
        per_snp_ancestry_background[name] = bm
    end

    # check correctness
    bm_tot = zeros(Int, nsamples, nsnps)
    for bm in values(per_snp_ancestry_background)
        bm_tot[bm] .+= 1
    end
    findmax(bm_tot)[1] == 1 || 
        error("sample $(findmax(bm_tot)[2][1]) SNP $(findmax(bm_tot)[2][2]) has >1 ancestry background, shouldn't happen")
    z = count(iszero, bm_tot)
    z > 0 && @warn("$z genotypes has no ancestry background identified")

    return per_snp_ancestry_background
end

function get_ancestry_specific_r2(
    genotype_file::String, # VCF (will import GT field)
    imputed_file::String, # VCF (import DS field, unless use_dosage=false, in which case we import GT)
    msp_file::String; # output of Rfmix2 (input to Rfmix2 MUST be genotype_file)
    ancestry_names::Vector{String} = get_ancestry_names(msp_file),
    summary_file::String = "", # must contain CHR/POS/REF/ALT/AF/isImputed columns
    use_dosage::Bool=true, # whether to import imputed data as dosage or genotypes
    maf_bins::Vector{Float64}=[0.0, 0.0005, 0.001, 0.005, 0.01, 0.02, 0.03, 0.04, 0.05, 0.1, 0.2, 0.3, 0.4, 0.5],
    )
    # import the genotype and imputed data matrices, after matching them
    Xgeno, Ximpt, mafs, samples, snps = import_Xtrue_Ximp(
        genotype_file, imputed_file, summary_file=summary_file, 
        use_dosage=use_dosage, 
    )

    # compute background ancestries on the phased genotypes
    n_ancestries = length(ancestry_names)
    ancestry_masks = create_LAI_mapping_matrices(
        samples, snps, msp_file, n_ancestries, ancestry_names=ancestry_names)

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

            # compute accuracy
            mask = has_this_ancestry[:, idx]
            truth = @view(Xgeno[:, idx][mask])
            imptd = @view(Ximpt[:, idx][mask])
            push!(R2, aggregate_R2(truth, imptd))
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

# include("/u/home/b/biona001/bge_analysis/concordance/scripts/utilities.jl")
# include("/u/home/b/biona001/bge_analysis/concordance/scripts/plots.jl")
# genotype_file = ["/u/home/b/biona001/project-loes/ForBen_genotypes_subset/QC_hg38_conformed_king/chr21.vcf.gz", "/u/home/b/biona001/project-loes/ForBen_genotypes_subset/QC_hg38_conformed_king/chr22.vcf.gz"]
# imputed_file = ["/u/home/b/biona001/project-loes/GLIMPSE2_toni/typedSNPs/chr21.sampleQC.snpQC.vcf.gz", "/u/home/b/biona001/project-loes/GLIMPSE2_toni/typedSNPs/chr22.sampleQC.snpQC.vcf.gz"]
# msp_file = ["/u/home/b/biona001/project-loes/ForBen_genotypes_subset/LAI/output_v2/chr21.msp.tsv", "/u/home/b/biona001/project-loes/ForBen_genotypes_subset/LAI/output_v2/chr22.msp.tsv"]
# summary_file = ["", ""]
# ancestry_names = get_ancestry_names(msp_file[1])
# use_dosage = true
# maf_bins = [0.0, 0.001, 0.005, 0.01, 0.02, 0.03, 0.04, 0.05, 0.1, 0.2, 0.3, 0.4, 0.5]
function get_ancestry_specific_r2(
    genotype_file::Vector{String}, # VCF files (will import GT field), e.g. one for each chr
    imputed_file::Vector{String}, # VCF files (import DS field, unless use_dosage=false, in which case we import GT)
    msp_file::Vector{String}; # output of Rfmix2 (input to Rfmix2 MUST be genotype_file)
    ancestry_names::Vector{String} = get_ancestry_names(msp_file[1]),
    summary_file::Vector{String} = ["" for _ in eachindex(genotype_file)], 
    use_dosage::Bool=true,
    maf_bins::Vector{Float64}=[0.0, 0.0005, 0.001, 0.005, 0.01, 0.02, 0.03, 0.04, 0.05, 0.1, 0.2, 0.3, 0.4, 0.5],
    )
    n_ancestries = length(ancestry_names)
    n_files = length(genotype_file)
    n_files == length(imputed_file) == length(msp_file) == 
        length(summary_file) || error("Input files should have same length")

    #
    # import everything (note: memory intensive if many chromosomes)
    # 
    ancestry_masks = []
    Xgenos = Matrix{Union{Missing, Float64}}[]
    Ximpts = Matrix{Union{Missing, Float64}}[]
    MAFs = Vector{Float64}[]
    shared_samples = Vector{String}[]
    shared_snps = Vector{SNP}[]
    for i in 1:n_files
        println("processing file $i / $n_files"); flush(stdout)
        # import the genotype and imputed data matrices, after matching them
        Xgeno, Ximpt, MAF, shared_sample, shared_snp = import_Xtrue_Ximp(
            genotype_file[i], imputed_file[i], summary_file=summary_file[i], 
            use_dosage=use_dosage, 
        )

        # compute background ancestries on the phased genotypes
        ancestry_mask = create_LAI_mapping_matrices(
            shared_sample, shared_snp, msp_file[i], n_ancestries, 
            ancestry_names=ancestry_names)

        push!(ancestry_masks, ancestry_mask)
        push!(Xgenos, Xgeno)
        push!(Ximpts, Ximpt)
        push!(MAFs, MAF)
        push!(shared_samples, shared_sample)
        push!(shared_snps, shared_snp)
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
    for ancestry in keys(ancestry_masks[1])
        df[!, ancestry] = zeros(size(df, 1)) # allocate R2 vector
    end

    #
    # aggregate LAI segments for each MAF bin
    #
    # pseudo code:
    # for EACH_MAF_BIN
    #     truth, imptd = Float64[], Float64[]
    #     for EACH_CHR 
    #         for EACH_ANCESTRY
    #             append_to_truth
    #             append_to_imptd
    #         end
    #     end
    #     R2 = aggregate_R2(truth, imptd)
    #     append_R2_to_df
    # end
    nsnps = zeros(length(maf_bins) - 1, length(ancestry_masks[1]))
    @showprogress for i in 2:length(maf_bins)
        bin_low, bin_high = maf_bins[i-1], maf_bins[i]
        truth = [Union{Missing, Float64}[] for _ in 1:length(keys(ancestry_masks[1]))]
        imptd = [Union{Missing, Float64}[] for _ in 1:length(keys(ancestry_masks[1]))]
        for j in 1:n_files
            # current Xtrue and Ximputed
            Xgeno, Ximpt, MAF, ancestry_mask = Xgenos[j], Ximpts[j], MAFs[j], ancestry_masks[j]
            idx = findall(x -> bin_low ≤ x < bin_high, MAF)
            length(idx) == 0 && continue

            # loop through different LAI
            for (k, has_this_ancestry) in enumerate(values(ancestry_mask))
                mask = has_this_ancestry[:, idx]
                append!(truth[k], @view(Xgeno[:, idx][mask]))
                append!(imptd[k], @view(Ximpt[:, idx][mask]))
            end
        end
        for (k, ancestry) in enumerate(keys(ancestry_masks[1]))
            println("Ancestry $ancestry has $(length(truth[k])) genotypes")
            df[i-1, ancestry] = aggregate_R2(truth[k], imptd[k])
            nsnps[i-1, k] += length(truth[k])
        end
    end

    # returns a dataframe containing aggregate R2 for different ancestries
    # the second return argument is the number of SNPs in each bin used to 
    # calculate the given R2
    return df, nsnps
end

# chr = 22
# genotype_file = "/u/home/b/biona001/project-loes/gpc_gsa/array/chr$chr.vcf.gz"
# imputed_file = "/u/home/b/biona001/project-loes/gpc_gsa/imputed/typedSNPs/chr22.vcf.gz"
# msp_file = "/u/home/b/biona001/project-loes/gpc_gsa/tbb/gpc_gsa_w2.chr22.typed.deconvolved.msp.tsv"
# summary_file = ""
# ancestry_names = get_ancestry_names(msp_file[1])
# summary_file = ""
# ancestry_names = get_ancestry_names(msp_file)
# maf_bins = collect(LinRange(0, 0.5, 100))
# df, nsnps = get_ancestry_specific_concordance(genotype_file, imputed_file, msp_file, maf_bins=maf_bins)
function get_ancestry_specific_concordance(
    genotype_file::String, # VCF (will import GT field)
    imputed_file::String, # VCF (import DS field, unless use_dosage=false, in which case we import GT)
    msp_file::String; # output of Rfmix2 (input to Rfmix2 MUST be genotype_file)
    ancestry_names::Vector{String} = get_ancestry_names(msp_file),
    summary_file::String = "", # must contain CHR/POS/REF/ALT/AF/isImputed columns
    maf_bins::Vector{Float64}=[0.0, 0.0005, 0.001, 0.005, 0.01, 0.02, 0.03, 0.04, 0.05, 0.1, 0.2, 0.3, 0.4, 0.5],
    )
    # import the genotype and imputed data matrices, after matching them
    Xgeno, Ximpt, mafs, samples, snps = import_Xtrue_Ximp(
        genotype_file, imputed_file, summary_file=summary_file, 
        use_dosage=false, 
    )

    # compute background ancestries on the phased genotypes
    n_ancestries = length(ancestry_names)
    ancestry_masks = create_LAI_mapping_matrices(
        samples, snps, msp_file, n_ancestries, ancestry_names=ancestry_names)

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
        # compute concordance
        nonref_concordances = Float64[]
        for i in 2:length(maf_bins)
            # find SNPs in current bin
            bin_low, bin_high = maf_bins[i-1], maf_bins[i]
            idx = findall(x -> bin_low ≤ x < bin_high, mafs)
            if length(idx) == 0
                @warn("maf_bin [$bin_low, $bin_high] has no SNP! Assining concordance of NaN")
                push!(nonref_concordances, NaN)
                continue
            end

            # compute accuracy
            mask = has_this_ancestry[:, idx]
            truth = @view(Xgeno[:, idx][mask])
            imptd = @view(Ximpt[:, idx][mask])
            contingency_table = compute_concordance(truth, imptd)
            push!(nonref_concordances, nonref_concordance(contingency_table))
            nsnps[i-1, counter] = length(idx) # record number of SNPs used to compute R2
        end
        df[!, ancestry] = nonref_concordances
        counter += 1
    end

    # returns a dataframe containing aggregate nonref_concordances for different ancestries
    # the second return argument is the number of SNPs in each bin used to 
    # calculate the given nonref_concordances
    return df, nsnps
end

# genotype_file = ["/u/home/b/biona001/project-loes/gpc_gsa/array/chr21.vcf.gz", "/u/home/b/biona001/project-loes/gpc_gsa/array/chr22.vcf.gz"]
# imputed_file = ["/u/home/b/biona001/project-loes/gpc_gsa/imputed/typedSNPs/chr21.vcf.gz", "/u/home/b/biona001/project-loes/gpc_gsa/imputed/typedSNPs/chr22.vcf.gz"]
# msp_file = ["/u/home/b/biona001/project-loes/gpc_gsa/tbb/gpc_gsa_w2.chr21.typed.deconvolved.msp.tsv", "/u/home/b/biona001/project-loes/gpc_gsa/tbb/gpc_gsa_w2.chr22.typed.deconvolved.msp.tsv"]
# genotype_file = ["/u/home/b/biona001/project-loes/ForBen_genotypes_subset/QC_hg38_conformed_king/chr21.vcf.gz", "/u/home/b/biona001/project-loes/ForBen_genotypes_subset/QC_hg38_conformed_king/chr22.vcf.gz"]
# imputed_file = ["/u/home/b/biona001/project-loes/GLIMPSE2_toni/typedSNPs/chr21.sampleQC.snpQC.vcf.gz", "/u/home/b/biona001/project-loes/GLIMPSE2_toni/typedSNPs/chr22.sampleQC.snpQC.vcf.gz"]
# msp_file = ["/u/home/b/biona001/project-loes/ForBen_genotypes_subset/LAI/output_v2/chr21.msp.tsv", "/u/home/b/biona001/project-loes/ForBen_genotypes_subset/LAI/output_v2/chr22.msp.tsv"]
# summary_file = ["", ""]
# ancestry_names = get_ancestry_names(msp_file[1])
# maf_bins = collect(LinRange(0, 0.5, 100))
# df, nsnps = get_ancestry_specific_concordance(genotype_file, imputed_file, msp_file, maf_bins=maf_bins)
function get_ancestry_specific_concordance(
    genotype_file::Vector{String}, # VCF files (will import GT field), e.g. one for each chr
    imputed_file::Vector{String}, # VCF files (import DS field, unless use_dosage=false, in which case we import GT)
    msp_file::Vector{String}; # output of Rfmix2 (input to Rfmix2 MUST be genotype_file)
    ancestry_names::Vector{String} = get_ancestry_names(msp_file[1]),
    summary_file::Vector{String} = ["" for _ in eachindex(genotype_file)], 
    maf_bins::Vector{Float64}=[0.0, 0.0005, 0.001, 0.005, 0.01, 0.02, 0.03, 0.04, 0.05, 0.1, 0.2, 0.3, 0.4, 0.5],
    )
    n_ancestries = length(ancestry_names)
    n_files = length(genotype_file)
    n_files == length(imputed_file) == length(msp_file) == 
        length(summary_file) || error("Input files should have same length")

    #
    # import everything (note: memory intensive if many chromosomes)
    # 
    ancestry_masks = []
    Xgenos = Matrix{Union{Missing, Float64}}[]
    Ximpts = Matrix{Union{Missing, Float64}}[]
    MAFs = Vector{Float64}[]
    shared_samples = Vector{String}[]
    shared_snps = Vector{SNP}[]
    for i in 1:n_files
        println("processing file $i / $n_files"); flush(stdout)
        # import the genotype and imputed data matrices, after matching them
        Xgeno, Ximpt, MAF, shared_sample, shared_snp = import_Xtrue_Ximp(
            genotype_file[i], imputed_file[i], summary_file=summary_file[i], 
            use_dosage=false, 
        )

        # compute background ancestries on the phased genotypes
        ancestry_mask = create_LAI_mapping_matrices(
            shared_sample, shared_snp, msp_file[i], n_ancestries, 
            ancestry_names=ancestry_names)

        push!(ancestry_masks, ancestry_mask)
        push!(Xgenos, Xgeno)
        push!(Ximpts, Ximpt)
        push!(MAFs, MAF)
        push!(shared_samples, shared_sample)
        push!(shared_snps, shared_snp)
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
    for ancestry in keys(ancestry_masks[1])
        df[!, ancestry] = zeros(size(df, 1)) # allocate R2 vector
    end

    #
    # aggregate LAI segments for each MAF bin
    #
    # pseudo code:
    # for EACH_MAF_BIN
    #     truth, imptd = Float64[], Float64[]
    #     for EACH_CHR 
    #         for EACH_ANCESTRY
    #             append_to_truth
    #             append_to_imptd
    #         end
    #     end
    #     concordance = compute_concordance(truth, imptd)
    #     append_concordance_to_df
    # end
    nsnps = zeros(length(maf_bins) - 1, length(ancestry_masks[1]))
    @showprogress for i in 2:length(maf_bins)
        bin_low, bin_high = maf_bins[i-1], maf_bins[i]
        truth = [Union{Missing, Float64}[] for _ in 1:length(keys(ancestry_masks[1]))]
        imptd = [Union{Missing, Float64}[] for _ in 1:length(keys(ancestry_masks[1]))]
        for j in 1:n_files
            # current Xtrue and Ximputed
            Xgeno, Ximpt, MAF, ancestry_mask = Xgenos[j], Ximpts[j], MAFs[j], ancestry_masks[j]
            idx = findall(x -> bin_low ≤ x < bin_high, MAF)
            length(idx) == 0 && continue

            # loop through different LAI
            for (k, has_this_ancestry) in enumerate(values(ancestry_mask))
                mask = has_this_ancestry[:, idx]
                append!(truth[k], @view(Xgeno[:, idx][mask]))
                append!(imptd[k], @view(Ximpt[:, idx][mask]))
            end
        end
        for (k, ancestry) in enumerate(keys(ancestry_masks[1]))
            println("Ancestry $ancestry has $(length(truth[k])) genotypes")

            # compute sensitivities/precision/non-ref-concordance for this ancestry
            contigency_table = compute_concordance(truth[k], imptd[k])
            df[i-1, ancestry] = nonref_concordance(contigency_table)
            nsnps[i-1, k] += length(truth[k])
            # sensitivities = sensitivity(contigency_table)
            # precisions = precision(contigency_table)
        end
    end

    # returns a dataframe containing nonref concordance for different ancestries
    # the second return argument is the number of SNPs in each bin used to 
    # calculate the given nonref_concordances
    return df, nsnps
end

"""
    compute_concordance(ximp::AbstractMatrix, xtrue::AbstractMatrix)

Computes sensitivity, precision, and non-reference concordances. Assume
we create the following contingency table:

                                        Imputed
                    |      ALT (1)       |       REF (0)       |
                    |------------------------------------------|
            ALT (1) | true positive  (A) | false negative (B)  |
Truth:              |------------------------------------------|
            REF (0) | false positive (C) | true negative (D)  |
                    |------------------------------------------|

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

+ Sensitivity = true positives / (true positives + false negaitves) = A / (A+B)
    (NOTE that the total number of ALT alleles in the array data is A+B)
+ Precision = true positives / (true positives + false positives) = A / (A+C)
+ Non-reference concordance = A / (A+B+C). To see why, consider the denominator
    to be the total number of ALT alleles called the *either* array/truth dataset
    or imputed dataset, i.e. (A+B+C). The numerator is where they agree, i.e. A.

# Returns
A vector of `ContingencyTable`, one for each SNP. Contingency tables support 
`sensitivity`, `precision`, and `nonref_concordance` functions, e.g. 
`sensitivity(x::ContingencyTable)`
"""
function compute_concordance(ximp::AbstractMatrix, xtrue::AbstractMatrix)
    size(ximp) == size(xtrue) || error("check dimension")
    n, p = size(ximp)
    contigency_tables = ContingencyTable[]
    for j in 1:p
        snp_j_contigency = compute_concordance(@view(ximp[:, j]), @view(xtrue[:, j]))
        push!(contigency_tables, snp_j_contigency)
    end
    return contigency_tables
end

function compute_concordance(ximp::AbstractVector, xtrue::AbstractVector)
    n = length(ximp)
    n == length(xtrue) || error("Expected ximp and xtrue to have the same length")
    A, B, C, D = 0, 0, 0, 0
    for i in 1:n
        ismissing(xtrue[i]) && continue
        ismissing(ximp[i]) && continue
        if xtrue[i] == 0 && ximp[i] == 0
            D += 2
        elseif xtrue[i] == 1 && ximp[i] == 0
            D += 1
            B += 1
        elseif xtrue[i] == 2 && ximp[i] == 0
            B += 2
        elseif xtrue[i] == 0 && ximp[i] == 1
            D += 1
            C += 1
        elseif xtrue[i] == 1 && ximp[i] == 1
            A += 1
            D += 1
        elseif xtrue[i] == 2 && ximp[i] == 1
            A += 1
            B += 1
        elseif xtrue[i] == 0 && ximp[i] == 2
            C += 2
        elseif xtrue[i] == 1 && ximp[i] == 2
            A += 1
            C += 1
        elseif xtrue[i] == 2 && ximp[i] == 2
            A += 2
        else
            error(
                "Expected xtrue[i,j] and ximpt[i,j] to take values " * 
                "0, 1, or 2. Did you call convert_gt on both?"
            )
        end
    end
    return ContingencyTable(A, B, C, D)
end

function compute_concordance(
    genotype_file::String, # VCF (will import GT field)
    imputed_file::String; # VCF (import DS field, unless use_dosage=false, in which case we import GT)
    summary_file::String = ""
    )

    # import the full genotype and imputed data matrices
    Xgeno, Ximpt, mafs, _, _ = import_Xtrue_Ximp(
        genotype_file, imputed_file, summary_file=summary_file, use_dosage=false,
    )

    # compute 2 by 2 contingency tables, one for each SNP
    contigency_tables = compute_concordance(Ximpt, Xgeno)

    # compute sensitivities/precision/non-ref-concordance
    sensitivities = sensitivity.(contigency_tables)
    precisions = precision.(contigency_tables)
    nonref_concordances = nonref_concordance.(contigency_tables)

    return sensitivities, precisions, nonref_concordances, mafs, contigency_tables
end

function compute_concordance(
    genotype_file::Vector{String}, # VCF (will import GT field)
    imputed_file::Vector{String}; # VCF (import DS field, unless use_dosage=false, in which case we import GT)
    summary_file::Vector{String} = ["" for _ in eachindex(genotype_file)]
    )
    contigency_tables, mafs = ContingencyTable[], Float64[]
    for (gtfile, impfile, summ_file) in zip(genotype_file, imputed_file, summary_file)
        # import the full genotype and imputed data matrices
        Xgeno, Ximpt, mf, _, _ = import_Xtrue_Ximp(
            gtfile, impfile, summary_file=summ_file, use_dosage=false,
        )

        # compute confordance for current file
        contigency_tables_subset = compute_concordance(Ximpt, Xgeno)
        append!(contigency_tables, contigency_tables_subset)
        append!(mafs, mf)
    end
    length(contigency_tables) == length(mafs) || error("length different!")

    # compute sensitivities/precision/non-ref-concordance
    sensitivities = sensitivity.(contigency_tables)
    precisions = precision.(contigency_tables)
    nonref_concordances = nonref_concordance.(contigency_tables)

    return sensitivities, precisions, nonref_concordances, mafs, contigency_tables
end
