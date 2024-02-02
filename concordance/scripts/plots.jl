using Plots
gr(fmt=:png)
include("utilities.jl")

function make_aggregate_R2_plot(
    aggregate_R2::Vector{T},
    plt_label::String,
    maf_bins::Vector{T},
    plt_title::String;
    ymin = 0.0, 
    ymax=1.0,
    markersize = ones(length(maf_bins) - 1),
    my_color = "mediumblue", # default color
    my_series_style = :solid, # default series style
    ) where T
    maf_bin_xaxis = [(maf_bins[i] + maf_bins[i-1])/2 for i in 2:length(maf_bins)]
    plt = plot(maf_bin_xaxis, aggregate_R2, markershape=:circle, 
        ylabel="Aggregate R2", label=plt_label,
        xlabel="MAF", ylim=(ymin, ymax), 
        # xaxis=:log, xlim=(0.00003, 0.6),
        # xticks=([0.0001, 0.001, 0.01, 0.1, 0.5], ["0.01", "0.1", "1", "10", "50"]),
        legend=:bottomright, title=plt_title,
        size=(400, 400), dpi=300, markerstrokewidth=0, w=3, 
        markersize=markersize, my_color=my_color, my_series_style=my_series_style
        )
    return plt
end

function make_aggregate_R2_plots(
    aggregate_R2s::Vector{Vector{T}},
    plt_labels::Vector{String},
    maf_bins::Vector{T},
    plt_title::String;
    ymin = 0.0, 
    ymax=1.0,
    markersize = ones(length(maf_bins) - 1),
    my_color = "mediumblue", # default color
    my_series_style = :solid, # default series style
    ) where T
    maf_bin_xaxis = [(maf_bins[i] + maf_bins[i-1])/2 for i in 2:length(maf_bins)]   

    # make init plot
    plt = make_aggregate_R2_plot(
        aggregate_R2s[1], plt_labels[1], maf_bins, plt_title,
        ymin=ymin, ymax=ymax, markersize=markersize[1],
        my_color=my_color[1], my_series_style=my_series_style[1]
    )

    # remaining plots
    for i in 2:length(aggregate_R2s)
        plot!(plt, maf_bin_xaxis, aggregate_R2s[i], label=plt_labels[i], 
            markershape=:circle, size=(400, 400), dpi=300, 
            markerstrokewidth=0, w=3, color=my_colors[i], ls=my_series_style[i]
            )
    end

    return plt
end

function make_concordance_plots(
    concordance::Vector{T}, # vector of non-ref concordance, one value for each SNP
    mafs::Vector{T}, # vector of MAF, one value for each SNP
    plt_label::String,
    maf_bins::Vector{T},
    plt_title::String;
    ymin = 0.0, 
    ymax = 1.0,
    markersize = ones(length(maf_bins) - 1),
    my_color = "mediumblue", # default color
    my_series_style = :solid, # default series style
    ) where T
    length(concordance) == length(mafs) || 
        error("expected length(concordance) == length(mafs)")
    maf_bin_xaxis = [(maf_bins[i] + maf_bins[i-1])/2 for i in 2:length(maf_bins)]

    # compute concordance by maf bins
    concordance_by_maf_bins = T[]
    for i in 2:length(maf_bins)
        idx = findall(x -> maf_bins[i-1] ≤ x ≤ maf_bins[i], mafs)
        println("MAF ∈ ($(maf_bins[i-1]), $(maf_bins[i])) have $(length(idx)) SNPs")
        cur_concordance = concordance[idx]
        non_nan_idx = findall(!isnan, cur_concordance)
        push!(concordance_by_maf_bins, mean(cur_concordance[non_nan_idx]))
    end

    plt = plot(maf_bin_xaxis, concordance_by_maf_bins, markershape=:circle, 
        ylabel="Non-ref conconcordance", label=plt_label,
        xlabel="MAF", ylim=(ymin, ymax), 
#         xaxis=:log, xlim=(0.0003, 0.6),
#         xticks=([0.001, 0.01, 0.1, 0.5], ["0.1", "1", "10", "50"]),
        legend=:bottomright, title=plt_title, markersize=markersize, 
        size=(400, 400), dpi=300, markerstrokewidth=0, w=3, 
        color=my_color, ls=my_series_style)

    return plt
end

function make_concordance_plots(
    concordances::Vector{Vector{T}},
    mafs::Vector{Vector{T}},
    plt_labels::Vector{String},
    maf_bins::Vector{T},
    plt_title::String;
    ymin = 0.0, 
    ymax = 1.0,
    markersizes = [ones(length(maf_bins) - 1) for _ in eachindex(concordances)],
    my_colors = palette(:default), # default series color (max length = 16)
    my_series_style = [:solid for _ in eachindex(concordances)] # default series style
    ) where T
    length(concordances) == length(mafs) == length(plt_labels) || 
        error("expected length(concordances) == length(mafs)")

    maf_bin_xaxis = [(maf_bins[i] + maf_bins[i-1])/2 for i in 2:length(maf_bins)]

    # make init plot
    plt = make_concordance_plots(
        concordances[1], mafs[1], plt_labels[1], maf_bins, plt_title,
        ymin=ymin, ymax=ymax, markersize=markersizes[1], 
        my_color=my_colors[1], my_series_style=my_series_style[1]
    )

    # remaining plots
    for j in 2:length(concordances)
        concordance = concordances[j]
        maf = mafs[j]

        concordance_by_maf_bins = T[]
        for i in 2:length(maf_bins)
            idx = findall(x -> maf_bins[i-1] ≤ x ≤ maf_bins[i], maf)
            cur_concordance = concordance[idx]
            non_nan_idx = findall(!isnan, cur_concordance)
            push!(concordance_by_maf_bins, mean(cur_concordance[non_nan_idx]))
        end

        plot!(plt, maf_bin_xaxis, concordance_by_maf_bins,
            label=plt_labels[j], markershape=:circle, 
            markersize=markersizes[j], size=(400, 400), dpi=300, 
            markerstrokewidth=0, w=3, color=my_colors[j], ls=my_series_style[j])
    end

    return plt
end

"""
    compute_bin_sizes(mafs, ngenotypes, maf_bins)

Helper function for computing bin sizes (intended to be used with
concordance plot by local ancestries)
"""
function compute_bin_sizes(
    mafs::AbstractVector, 
    ngenotypes::AbstractVector,
    maf_bins::AbstractVector;
    offset = 8 # scaling factor to make sizes look reasonable
    )
    bin_sizes = Int[]

    for i in 2:length(maf_bins)
        # first count number of SNPs in each maf bin
        idx = findall(x -> maf_bins[i-1] ≤ x ≤ maf_bins[i], mafs)
        nsnps = length(idx)

        # scale bin sizes by number of genotypes with given local ancestry background
        push!(bin_sizes, nsnps * sum(ngenotypes[idx]))
    end

    bin_sizes_normalized = log.(bin_sizes) .- offset
    return bin_sizes_normalized
end
