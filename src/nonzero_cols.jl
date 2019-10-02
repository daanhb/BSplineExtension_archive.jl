
using ..CompactInfiniteVectors, InfiniteVectors
using GridArrays.ModCartesianIndicesBase: ModCartesianIndices
using FrameFun.BasisFunctions: Dictionary, Measure, AbstractGrid, grid, support
using GridArrays: mask


export nonzero_cols

nonzero_cols(dict::Dictionary, measure::Measure) =
    nonzero_cols(dict, grid(measure))
function nonzero_cols(dict::Dictionary, grid::AbstractGrid)
    q = div.(size(supergrid(grid)),size(dict))
    r = rem.(size(supergrid(grid)),size(dict))
    @assert all(r.==0)
    nonzero_cols(compactinfinitevector(dict, grid), q, mask(grid))
end

nonzero_cols(b::CompactInfiniteVector, m::Tuple{Int}, gridmask::BitArray{1}) =
    (s=InfiniteVectors.support(b);nonzero_cols(CartesianIndex(first(s)):CartesianIndex(last(s)), m, gridmask))

nonzero_cols(b::NTuple{N,CompactInfiniteVector}, m::NTuple{N,Int}, gridmask::BitArray{N}) where N =
    nonzero_cols(support(b), m, gridmask)

function nonzero_cols(b_support::CartesianIndices{N}, m::NTuple{N,Int}, gridmask::BitArray{N}) where N
    gridsize = size(gridmask)
    dictsize = div.(gridsize,m)
    A = falses(dictsize)
    for k in CartesianIndices(dictsize)
        l = (k.I .- 1).*m
        bk_support = b_support .+ CartesianIndex(l)
        support_in = false
        support_out = false
        for i in ModCartesianIndices(gridsize, first(bk_support), last(bk_support))
            if gridmask[i]
                support_in = true
            else
                support_out = true
            end
            if support_in && support_out
                A[k] = true
                break;
            end
        end
    end
    findall(A)
end




# # TODO
# # nonzero_cols(dict::TensorProductDict, gamma::GridArrays.ProductGrid, domain::Domains.ProductDomain) =


# export nonzero_coefficients
# """
#     nonzero_coefficients(dict::Dictionary1d, x::Real)
#
# Return the coefficient indices (`UnitRange`) of the `dict` elements that evaluate non zero in `x`
# """
# function nonzero_coefficients(dict::Dictionary1d, x::Real)
#     # The init_index is the starting index of all spline elements that overlap with x
#     tuple(1:length(dict))
# end
# #
# # # Limits of the indices of the coefficients of B that overlap with x.
# # # This is a tupple of (number of elements in tuple depends is equal to dimension) of two element tuples.
# # nonzero_coefficients(B::Dictionary, x::Real) =
# #     tuple(nonzero_coefficients(B, x))
#
# nonzero_coefficients(dict::TensorProductDict{N}, x) where N =
#     ntuple(k->nonzero_coefficients(elements(dict)[k],x[k])[1] ,Val(N))
#
# # Cartesian index limits of the coefficients of B that overlap with x.
# function coefficient_cartesian_index_limits_of_overlapping_elementst(B::Dictionary, x)
#     range = nonzero_coefficients(B, x)
#     CartesianIndex(range[1][1]), CartesianIndex(range[1][end])
# end
#
# function coefficient_cartesian_index_limits_of_overlapping_elementst(B::Dictionary, x::SVector{N}) where N
#     ranges = nonzero_coefficients(B, x)
#     CartesianIndex(ntuple(k->ranges[k][1] ,Val(N))), CartesianIndex(ntuple(k->ranges[k][end] ,Val(N)))
# end
#
# # Range of coefficient indices of B that overlap with the point x.
# coefficient_index_range_of_overlapping_elements(B::Dictionary, x) =
#     ModCartesianIndices(size(B), coefficient_cartesian_index_limits_of_overlapping_elementst(B, x)...)
#
# coefficient_index_mask_of_overlapping_elements(d::Dictionary, g::AbstractArray) =
#     coefficient_index_mask_of_overlapping_elements!(BitArray(undef, size(d)), d, g)
#
# function coefficient_index_mask_of_overlapping_elements!(mask::BitArray, B::Dictionary, g::AbstractArray)
#     fill!(mask, 0)
#     for x in g, i in coefficient_index_range_of_overlapping_elements(B, x)
#         mask[i] = 1
#     end
#     mask
# end
#
# coefficient_indices_of_overlapping_elements(dict::Dictionary, boundary::AbstractArray) =
#     findall(coefficient_index_mask_of_overlapping_elements(dict, boundary))
