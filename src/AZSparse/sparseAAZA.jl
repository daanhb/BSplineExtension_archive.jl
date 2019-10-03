

using ..BSplineExtensionSolvers: nonzero_cols, nonzero_rows
using InfiniteVectors, GridArrays, FrameFun.BasisFunctions, CardinalBSplines, CompactTranslatesDict, SparseArrays, FrameFun
using FrameFun.BasisFunctions: DiscreteMeasure, support
using GridArrays.ModCartesianIndicesBase: ModCartesianIndices
using ..CompactInfiniteVectors

function sparseAAZAmatrix(dict1::Dictionary, dict2::Dictionary, grid::AbstractGrid)
    ImZA = sparsemixedgramcomplement(dict2, dict1, discretemeasure(grid))
    col_indices = findall(reshape(nonzero_rows(ImZA),size(dict1)))
    RAE = sparseRAE(dict1, grid, col_indices)

    RAE*ImZA[LinearIndices(size(dict1))[col_indices],:]
end

sparsemixedgramcomplement(dict1::Dictionary, dict2::Dictionary, measure::DiscreteMeasure) =
    sparsemixedgramcomplement(nonzero_cols(dict2, measure), compactinfinitevector(dict1, grid(measure)),
        compactinfinitevector(dict2, grid(measure)), _oversamplingfactor(dict1,grid(measure)), grid(measure))

_oversamplingfactor(dict::ExtensionFrame, grid::AbstractGrid) =
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
        m::NTuple{N,Int}, dff::CartesianIndices{N}, grid::AbstractGrid) where N
    L = length(dff)

    R = Matrix{T}(undef, L, length(indices))
    zerocartesian = CartesianIndex(ntuple(k->0,Val(N)))

    bb = Array(OuterProductArray(map(subvector, b)...))
    b̃b̃ = Array(OuterProductArray(map(subvector, b̃)...))

    gridmask = mask(grid)
    gridsize = size(supergrid(grid))

    for (i,k) in enumerate(indices)
        for (j,d) in enumerate(dff)
            R[j,i] = - _innerproduct(bb, b̃b̃, b_support, b̃_support, m, k, d, gridmask, gridsize)
            if d==zerocartesian
                R[j,i] += 1
            end
        end
    end
    R
end

function sparsemixedgramcomplement(indices::AbstractArray{CartesianIndex{N}},
        b̃::NTuple{N,CompactInfiniteVector}, b::NTuple{N,CompactInfiniteVector}, m::NTuple{N,Int}, grid::AbstractGrid) where N
    T = promote_type(map(eltype,b)..., map(eltype, b̃)...)
    b_support = support(b)
    b̃_support = support(b̃)
    dff = difference_indices(b,b̃,m)

    b_supportlength = length(dff)

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
        ls = ModCartesianIndices(basissize,first(kdff),last(kdff))

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
        for j in 1:b_supportlength
            colix[j] = j
        end
        sort!(colix ,1,b_supportlength, InsertionSort,Base.Perm(Base.Order.ForwardOrdering(),rowvalscol))
        # sortperm!(colix,rowvalscol)
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
    grid::AbstractGrid, indices::AbstractVector{CartesianIndex{N}},
    dictsize::NTuple{N,Int}) where N
    gridsize = size(supergrid(grid))

    B = Array(OuterProductArray(map(subvector, b)...))[:]
    b_support = support(b)
    b_supportlength = length(b_support)

    m = div.(gridsize,dictsize)
    L = LinearIndices(gridsize)

    nnz = b_supportlength*length(indices)

    colptr = Vector{Int}(undef, length(indices)+1)
    nzvals = Vector{eltype(B)}(undef, nnz)
    rowvals = Vector{Int}(undef, nnz)

    rowvalscol = Vector{Int}(undef, b_supportlength)
    nzvalscol = Vector{eltype(B)}(undef, b_supportlength)
    colix = Vector{Int}(undef, b_supportlength)

    gridmask = mask(grid)
    newindices = cumsum(gridmask[:])

    colptr[1] = 1
    nzvalindex = 1
    for (i,k) in enumerate(indices)
        # support of element with index k
        bk_support = CartesianIndex(m.*(k.I.-1)) .+ b_support
        bk_support_indices = ModCartesianIndices(gridsize, first(bk_support), last(bk_support))
        colptrcol = 0
        for (j,l) in enumerate(bk_support_indices)
            if gridmask[l]
                colptrcol += 1
                rowvalscol[colptrcol] = newindices[L[l]]
                nzvalscol[colptrcol] = B[j]
            end
        end
        for j in colptrcol+1:b_supportlength
            rowvalscol[j] = length(L)+1
        end

        for j in 1:b_supportlength
            colix[j] = j
        end
        sort!(colix ,1,b_supportlength, InsertionSort,Base.Perm(Base.Order.ForwardOrdering(),rowvalscol))
        # sortperm!(colix, rowvalscol)
        for j in 1:colptrcol
            rowvals[nzvalindex] = rowvalscol[colix[j]]
            nzvals[nzvalindex] = nzvalscol[colix[j]]
            nzvalindex += 1
        end
        colptr[i+1] = colptr[i] + colptrcol
    end

    SparseMatrixCSC(length(grid),length(indices),colptr,resize!(rowvals,nzvalindex-1),resize!(nzvals,nzvalindex-1))
end



sparseRAE(dict::Dictionary, grid::AbstractGrid, indices) =
    sparseRAE(compactinfinitevector(dict, grid), grid, indices, size(dict))

@inline function _innerproduct(bb::Array{T,N}, b̃b̃::Array{T,N}, b_support::CartesianIndices{N},
    b̃_support::CartesianIndices{N},
    m::NTuple{N,Int}, k::CartesianIndex{N}, coef_dff::CartesianIndex{N}, gridmask::BitArray{N}, gridsize::NTuple{N,Int}) where {N,T}

    grid_dff = CartesianIndex(m .* coef_dff.I)
    overlapping_support = overlapping(first(b̃_support) + grid_dff, last(b̃_support) + grid_dff,first(b_support),last(b_support))

    ktranslate = CartesianIndex((k.I .-1 ) .* m)

    overlapping_k_support = overlapping_support .+ ktranslate
    # prepare to put overlapping_k_support as index in maskedgrid
    I = GridArrays.ModCartesianIndicesBase.ModCartesianIndices(gridsize, first(overlapping_k_support), last(overlapping_k_support))

    r = T(0)
    @inbounds for (i,imod) in zip(overlapping_k_support,I)
        if gridmask[imod]
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
    overlapping(l1,r1,l2,r2)
end

function overlapping(l1,r1,l2,r2)
    if l1 <= l2 <= r1
        return l2:r1
    elseif l2 <= l1 <= r2
        return l1:r2
    else
        return r1:l1
    end

end
