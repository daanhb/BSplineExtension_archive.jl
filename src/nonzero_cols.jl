export nonzero_cols
nonzero_cols(dict::Dictionary, gamma::AbstractGrid, domain::Domain) =
    coefficient_indices_of_overlapping_elements(dict, gamma[findall(boundary_mask(gamma, domain, true))])

nonzero_cols(dict::TensorProductDict{N}, gamma::AbstractGrid, domain::Domain) where N =
    coefficient_indices_of_overlapping_elements(dict, gamma[findall(boundary_mask(gamma, domain, ntuple(k->true,Val(N))))])

# # TODO
# # nonzero_cols(dict::TensorProductDict, gamma::GridArrays.ProductGrid, domain::Domains.ProductDomain) =


export nonzero_coefficients
"""
    nonzero_coefficients(dict::Dictionary1d, x::Real)

Return the coefficient indices (`UnitRange`) of the `dict` elements that evaluate non zero in `x`
"""
function nonzero_coefficients(dict::Dictionary1d, x::Real)
    # The init_index is the starting index of all spline elements that overlap with x
    tuple(1:length(dict))
end
#
# # Limits of the indices of the coefficients of B that overlap with x.
# # This is a tupple of (number of elements in tuple depends is equal to dimension) of two element tuples.
# nonzero_coefficients(B::Dictionary, x::Real) =
#     tuple(nonzero_coefficients(B, x))

nonzero_coefficients(dict::TensorProductDict, x) =
    tuple(map((dicti,xi)->nonzero_coefficients(dicti,xi)[1], elements(dict), x)...)

# Cartesian index limits of the coefficients of B that overlap with x.
function coefficient_cartesian_index_limits_of_overlapping_elementst(B::Dictionary, x)
    range = nonzero_coefficients(B, x)
    CartesianIndex(range[1][1]), CartesianIndex(range[1][end])
end

function coefficient_cartesian_index_limits_of_overlapping_elementst(B::Dictionary, x::SVector{N}) where N
    ranges = nonzero_coefficients(B, x)
    CartesianIndex(ntuple(k->ranges[k][1] ,Val(N))), CartesianIndex(ntuple(k->ranges[k][end] ,Val(N)))
end

# Range of coefficient indices of B that overlap with the point x.
coefficient_index_range_of_overlapping_elements(B::Dictionary, x) =
    ModCartesianIndices(size(B), coefficient_cartesian_index_limits_of_overlapping_elementst(B, x)...)

coefficient_index_mask_of_overlapping_elements(d::Dictionary, g::AbstractArray) =
    coefficient_index_mask_of_overlapping_elements!(BitArray(undef, size(d)), d, g)

function coefficient_index_mask_of_overlapping_elements!(mask::BitArray, B::Dictionary, g::AbstractArray)
    fill!(mask, 0)
    for x in g, i in coefficient_index_range_of_overlapping_elements(B, x)
        mask[i] = 1
    end
    mask
end

coefficient_indices_of_overlapping_elements(dict::Dictionary, boundary::AbstractArray) =
    findall(coefficient_index_mask_of_overlapping_elements(dict, boundary))
