# Write fast and sparse I-Z'A creator
# Write sparse EB struct
# Write sparse multiplication
# make number of nonzerocolumns smaller
# check correctness of difference_indices

# Write fast and sparse I-Z'A creator
# Write something that returns rowvals, colvals, nzvals for R*M⊗...⊗M
# Write fast multiplier for the sparse structs above
# write nonzero_cols if domain is a grid

module AZSparse
using ..BSplineExtensionSolvers: nonzero_cols


using InfiniteVectors, GridArrays, FrameFun.BasisFunctions, CardinalBSplines, CompactTranslatesDict, SparseArrays, FrameFun
using InfiniteVectors: sublength, _firstindex, _lastindex
using GridArrays: MaskedGrid, AbstractSubGrid
using FrameFun.BasisFunctions: unsafe_matrix, DiscreteMeasure
using CompactTranslatesDict: PeriodicEquispacedTranslates
using GridArrays.ModCartesianIndicesBase: ModCartesianIndices

import FrameFun.BasisFunctions: support
import CardinalBSplines: bsplinesignal

sparsemixedgramcomplement(dict1::Dictionary, dict2::Dictionary, measure::DiscreteMeasure;options...) =
    sparsemixedgramcomplement(nonzero_cols(dict2, measure), bsplinesignal(dict1, measure),
        bsplinesignal(dict2, measure), _oversamplingfactor(dict1,measure), FrameFun.grid(measure))

_oversamplingfactor(dict::Dictionary, measure::DiscreteMeasure) =
    _oversamplingfactor(dict, FrameFun.grid(measure))
_oversamplingfactor(dict::ExtensionFrame, grid::AbstractSubGrid) =
    _oversamplingfactor(basis(dict), supergrid(grid))
_oversamplingfactor(dict::TensorProductDict, grid::ProductGrid) =
    map(_oversamplingfactor, elements(dict), elements(grid))
function _oversamplingfactor(dict::Dictionary, grid::AbstractIntervalGrid)
    q, r = divrem(length(grid), length(dict))
    @assert r == 0
    q
end

function difference_indices(b::NTuple{N,CompactInfiniteVector}, b̃::NTuple{N,CompactInfiniteVector}, m::NTuple{N,Int}) where N
    b_support = support(b)
    b̃_support = support(b̃)

    dff1 = first(b_support)-last(b̃_support)
    dff2 = last(b_support)-first(b̃_support)

    CartesianIndex(map((x,y)->trunc(Int,x/y), dff1.I, m)):CartesianIndex(map((x,y)->trunc(Int,x/y), dff2.I,m))
end

function sparsemixedgramcomplement_nzband(indices::AbstractArray{CartesianIndex{N}}, T,
        b̃::NTuple{N,CompactInfiniteVector}, b̃_support::CartesianIndices{N},
        b::NTuple{N,CompactInfiniteVector}, b_support::CartesianIndices{N},
        m::NTuple{N,Int}, dff::CartesianIndices{N}, grid::AbstractSubGrid) where N
    L = length(dff)

    R = Matrix{T}(undef, L, length(indices))
    zerocartesian = CartesianIndex(ntuple(k->0,Val(N)))
    for (i,k) in enumerate(indices)
        for (j,d) in enumerate(dff)
            R[j,i] = (d==zerocartesian ? T(1) : T(0)) - _innerproduct(b, b_support, b̃, b̃_support, m, k, d, grid)
        end
    end
    R
end

function sparsemixedgramcomplement(indices::AbstractArray{CartesianIndex{N}},
        b̃::NTuple{N,CompactInfiniteVector}, b::NTuple{N,CompactInfiniteVector}, m::NTuple{N,Int}, grid::AbstractSubGrid) where N
    T = promote_type(map(eltype,b)..., map(eltype, b̃)...)
    b_support = support(b)
    b̃_support = support(b̃)
    dff = difference_indices(b,b̃,m)

    R = sparsemixedgramcomplement_nzband(indices, T, b̃, b̃_support, b, b_support, m, dff, grid)

    nonzeromaskR = .!(abs.(R).+1 .≈ 1)
    nnz = count(nonzeromaskR)
    basissize = div.(size(supergrid(grid)),m)
    n = length(indices)
    m = prod(basissize)
    nzvals = Vector{T}(undef, nnz)
    colptr = Vector{Int}(undef, n+1)
    rowvals = Vector{Int}(undef, nnz)
    L = LinearIndices(basissize)
    # periodic = ntuple(k->true,Val(N))
    valn = 1
    colptr[1] = 1
    for (i,k) in enumerate(indices)
        nnzcol = 0
        # ls = ModCartesianIndices(basissize, k+first(dff),k+last(dff),periodic)

        kdff = k .+ dff
        ls = GridArrays.ModCartesianIndicesBase.ModCartesianIndices(basissize,first(kdff),last(kdff))
        @assert length(ls) == length(dff)
        for (j,l) in enumerate(ls)
            if nonzeromaskR[j,i]
                nzvals[valn] = R[j,i]
                rowvals[valn] = L[l]
                nnzcol += 1
                valn += 1
            end
        end
        colptr[i+1] = colptr[i] + nnzcol
    end
    SparseMatrixCSC(m,n,colptr,rowvals,nzvals)
end



support(vecs::NTuple{N,CompactInfiniteVector}) where N =
    CartesianIndex(map(_firstindex, vecs)):CartesianIndex(map(_lastindex, vecs))

product(b::NTuple{N,CompactInfiniteVector}, i::CartesianIndex{N}) where {N} =
    prod(map(getindex, b, i.I))
bsplinesignal(dict::Dictionary, measure::DiscreteMeasure) =
    bsplinesignal(dict, FrameFun.grid(measure))
bsplinesignal(dict::ExtensionFrame, grid::AbstractSubGrid) = bsplinesignal(basis(dict), supergrid(grid))
bsplinesignal(dict::TensorProductDict, grid::ProductGrid) = map(bsplinesignal, elements(dict), elements(grid))
function bsplinesignal(dict::PeriodicEquispacedTranslates{T}, grid::AbstractEquispacedGrid) where {T}
    A = evaluation_operator(dict, grid)
    @assert A isa VerticalBandedOperator
    a = convert(Vector{T}, unsafe_matrix(A).array)
    f = 1
    l = length(a)
    for i in 1:length(a)
        if !(a[i] + 1 ≈ 1)
            break
        end
        f += 1
    end
    for j in length(a):-1:1
        if !(a[j] + 1 ≈ 1)
            break
        end
        l -= 1
    end
    truncatedarray = a[f:l]
    os = mod(unsafe_matrix(A).offset+f, length(grid))
    if !(0 <= os + length(truncatedarray) <= length(grid))
        os -= length(grid)
    end
    CompactInfiniteVector(truncatedarray, os)
end

function _innerproduct(b::NTuple{N,CompactInfiniteVector{T}}, b_support::CartesianIndices{N},
    b̃::NTuple{N,CompactInfiniteVector{T}}, b̃_support::CartesianIndices{N},
    m::NTuple{N,Int}, k::CartesianIndex{N}, d::CartesianIndex{N}, maskedgrid::MaskedGrid) where {N,T}

    coef_dff = d
    grid_dff = CartesianIndex(m .* coef_dff.I)
    b̃dff_support = b̃_support .+ grid_dff
    overlapping_support = overlapping(b̃dff_support,b_support)
    @assert length(overlapping_support) > 0

    ktranslate = CartesianIndex((k.I .-1 ) .* m)

    overlapping_k_support = overlapping_support .+ ktranslate
    # prepare to put overlapping_k_support as index in maskedgrid
    I = GridArrays.ModCartesianIndicesBase.ModCartesianIndices(size(supergrid(maskedgrid)), first(overlapping_k_support), last(overlapping_k_support))

    r = T(0)
    for (i,imod) in zip(overlapping_k_support,I)
        if maskedgrid.mask[imod]
            r += product(b, i - ktranslate)*product(b̃, i-ktranslate-grid_dff)
        end
    end
    r
end

function overlapping(a::CartesianIndices, b::CartesianIndices)
    l1,r1 = first(a), last(a)
    l2,r2 = first(b), last(b)

    if l1 <= l2 <= r1
        return l2:r1
    elseif l2 <= l1 <= r2
        return l1:r2
    else
        return r1:l1
    end

end

##########################################
# using StaticArrays, BSplineExtension, FrameFun, LinearAlgebra
# D = .4*disk() + SVector(.5,.5)
# Pbasis = NdCDBSplinePlatform((3,3))
# P = ExtensionFramePlatform(Pbasis, D)
# N = 20,20
# dict1 = dictionary(P,N)
# g = oversampling_grid(P,N)
# γ = supergrid(g)
# μ = discretemeasure(g)
# dict2 = dualdictionary(P,N,μ)
# L = LinearIndices(size(dict1))
# C = CartesianIndices(size(dict1))
#
# G = Matrix(mixedgramoperator(dict2,dict1,μ))
# indices = nonzero_cols(dict1, μ)
# GG = copy(G);GG[:,L[indices]] .= 0;
# sum(abs.(GG))≈count(sum(abs.(GG),dims=1) .≈1 )
#
# IG = I-G
# IGE = IG[:,L[indices]]
#
# b = bsplinesignal(dict1,μ)
# b̃ = bsplinesignal(dict2,μ)
#
# dff = difference_indices(b,b̃,(2,2))
#
# R = Matrix{Float64}(undef, length(dff), length(indices))
#
# M = GridArrays.ModCartesianIndicesBase.ModCartesianIndices(size(basis(dict1)),first(dff), last(dff))
# M1 = GridArrays.ModCartesianIndicesBase.ModCartesianIndices(size(basis(dict1)),CartesianIndex(1,1),CartesianIndex(size(basis(dict1))))
# for (i,k) in enumerate(indices)
#     kdff = k .+ dff
#     M = GridArrays.ModCartesianIndicesBase.ModCartesianIndices(size(basis(dict1)),first(kdff),last(kdff))
#     for (j,l) in enumerate(M)
#         R[j,i] = IG[L[M1[l]],L[k]]
#     end
# end
# sum(abs.(R),dims=1) ≈ sum(abs.(IGE),dims=1)
# count(sum(abs.(R),dims=2) .+ 1 .≈ 1)
#
# RR = sparsemixedgramcomplement_nzband(indices, Float64, b̃, support(b̃), b, support(b), (2,2), dff, g)
# RR≈R
# s = sparsemixedgramcomplement(dict2,dict1,μ)
# sparse(IGE)≈s



end
