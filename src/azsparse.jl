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
using ..BSplineExtensionSolvers: nonzero_cols, nonzero_rows


using InfiniteVectors, GridArrays, FrameFun.BasisFunctions, CardinalBSplines, CompactTranslatesDict, SparseArrays, FrameFun
using InfiniteVectors: sublength, _firstindex, _lastindex
using GridArrays: MaskedGrid, AbstractSubGrid
using FrameFun.BasisFunctions: unsafe_matrix, DiscreteMeasure
using CompactTranslatesDict: PeriodicEquispacedTranslates
using GridArrays.ModCartesianIndicesBase: ModCartesianIndices

import FrameFun.BasisFunctions: support
import CardinalBSplines: bsplinesignal

function sparseAAZAmatrix(dict1, dict2, measure)
    ImZA = sparsemixedgramcomplement(dict2, dict1, measure)
    col_indices = findall(reshape(nonzero_rows(ImZA),size(dict1)))
    RAE = sparseRAE(dict1, measure, col_indices)

    RAE*ImZA[LinearIndices(size(dict1))[col_indices],:]
end

sparsemixedgramcomplement(dict1::Dictionary, dict2::Dictionary, measure::DiscreteMeasure) =
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

    bb = Array(OuterProductArray(map(subvector, b)...))
    b̃b̃ = Array(OuterProductArray(map(subvector, b̃)...))

    for (i,k) in enumerate(indices)
        for (j,d) in enumerate(dff)
            R[j,i] = (d==zerocartesian ? T(1) : T(0)) - _innerproduct(bb, b̃b̃, b_support, b̃_support, m, k, d, grid)
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

    rowvalscol = Vector{Int}(undef, length(dff))
    nzvalscol = Vector{T}(undef, length(dff))
    colix = Vector{Int}(undef, length(dff))

    L = LinearIndices(basissize)
    # periodic = ntuple(k->true,Val(N))
    valn = 1
    colptr[1] = 1
    @inbounds for (i,k) in enumerate(indices)
        nnzcol = 0
        kdff = k .+ dff
        ls = GridArrays.ModCartesianIndicesBase.ModCartesianIndices(basissize,first(kdff),last(kdff))

        for (j,l) in enumerate(ls)
            if nonzeromaskR[j,i]
                nnzcol += 1
                nzvalscol[nnzcol] = R[j,i]
                rowvalscol[nnzcol] = L[l]
            end
        end
        for j in nnzcol+1:length(dff)
            rowvalscol[j] = length(L)+1
        end
        sortperm!(colix, rowvalscol)
        for j in 1:nnzcol
            nzvals[valn] = nzvalscol[colix[j]]
            rowvals[valn] = rowvalscol[colix[j]]
            valn += 1
        end

        colptr[i+1] = colptr[i] + nnzcol
    end
    SparseMatrixCSC(m,n,colptr,rowvals,nzvals)
end

function sparseRAE(b::NTuple{N,CompactInfiniteVector},
    grid::MaskedGrid, indices::AbstractVector{CartesianIndex{N}},
    dictsize::NTuple{N,Int}) where N
    gridsize = size(supergrid(grid))

    B = Array(OuterProductArray(map(subvector, b)...))[:]
    b_support = support(b)
    b_supportlength = length(b_support)

    m = div.(gridsize,dictsize)
    L = LinearIndices(gridsize)

    nnz = b_supportlength*length(indices)

    colptr = collect(1:b_supportlength:nnz+1)
    nzvals = Vector{eltype(B)}(undef, nnz)
    rowvals = Vector{Int}(undef, nnz)

    rowvalscol = Vector{Int}(undef, b_supportlength)
    nzvalscol = Vector{Int}(undef, b_supportlength)
    colix = Vector{Int}(undef, b_supportlength)

    for (i,k) in enumerate(indices)
        # support of element with index k
        bk_support = CartesianIndex(m.*(k.I.-1)) .+ b_support
        bk_support_indices = ModCartesianIndices(gridsize, first(bk_support), last(bk_support))

        for (j,l) in enumerate(bk_support_indices)
            rowvalscol[j] = L[l]
        end
        sortperm!(colix, rowvalscol)
        j = (i-1)*b_supportlength + 1
        for (colixi) in (colix)
            rowvals[j] = rowvalscol[colixi]
            nzvals[j] = B[colixi]
            j += 1
        end
    end

    A = SparseMatrixCSC(prod(gridsize),length(indices),colptr,rowvals,nzvals)
    A[L[subindices(grid)],:]
end
sparseRAE(dict::Dictionary, measure::DiscreteMeasure, indices) =
    sparseRAE(bsplinesignal(dict, measure), grid(measure), indices, size(dict))

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

function _innerproduct(bb::Array{T,N}, b̃b̃::Array{T,N}, b_support::CartesianIndices{N},
    b̃_support::CartesianIndices{N},
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
    @inbounds for (i,imod) in zip(overlapping_k_support,I)
        if maskedgrid.mask[imod]
            k = i - ktranslate
            if k∈b_support
                b1 = bb[(k -first(b_support)+CartesianIndex{N}(1))]
                k = i-ktranslate-grid_dff
                if k∈b̃_support
                    b2 = b̃b̃[(k -first(b̃_support)+CartesianIndex{N}(1))]
                    r += b1*b2
                end
            end
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
end
