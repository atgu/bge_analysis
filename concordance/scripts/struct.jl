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
The `ContingencyTable` struct defines a 2 by 2 contigency table between an 
imputed SNP and the ground truth genotype. The fields A/B/C/D defines true
positives, false negatives, false positives, and true negatives. See the
function `compute_concordance` for how A/B/C/D are computed. 
"""
struct ContingencyTable
    A::Int # true positives
    B::Int # false negatives
    C::Int # false positives
    D::Int # true negatives
end
sensitivity(x::ContingencyTable) = x.A / (x.A + x.B)
precision(x::ContingencyTable) = x.A / (x.A + x.C)
nonref_concordance(x::ContingencyTable) = x.A / (x.A + x.B + x.C)
