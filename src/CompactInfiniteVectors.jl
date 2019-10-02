module CompactInfiniteVectors
using GridArrays: supergrid, AbstractEquispacedGrid, AbstractGrid, ProductGrid, AbstractSubGrid
using FrameFun.BasisFunctions: evaluation_operator, VerticalBandedOperator, unsafe_matrix, elements, TensorProductDict
using FrameFun: basis, ExtensionFrame
using InfiniteVectors: CompactInfiniteVector, sublength, _firstindex, _lastindex
using CompactTranslatesDict: PeriodicEquispacedTranslates

import FrameFun.BasisFunctions: support

support(vecs::NTuple{N,CompactInfiniteVector}) where N =
    CartesianIndex(map(_firstindex, vecs)):CartesianIndex(map(_lastindex, vecs))

product(b::NTuple{N,CompactInfiniteVector}, i::CartesianIndex{N}) where {N} =
    prod(map(getindex, b, i.I))

export compactinfinitevector
compactinfinitevector(dict::ExtensionFrame, grid::AbstractSubGrid) = compactinfinitevector(basis(dict), supergrid(grid))
compactinfinitevector(dict::TensorProductDict, grid::ProductGrid) = map(compactinfinitevector, elements(dict), elements(grid))
function compactinfinitevector(dict::PeriodicEquispacedTranslates{T}, grid::AbstractEquispacedGrid) where {T}
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
end
