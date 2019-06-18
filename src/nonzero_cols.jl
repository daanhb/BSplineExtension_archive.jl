
nonzero_cols(platform::FrameFun.Platform, param; dict=dictionary(platform, param), options...) =
    nonzero_cols(dict, sampling_grid(platform.basisplatform, param; dict=basis(dict), options...))

nonzero_cols(dict::ExtensionFrame, gamma::AbstractGrid) =
    nonzero_cols(basis(dict), gamma, support(dict))

nonzero_cols(dict::Dictionary, gamma::AbstractGrid, domain::Domain) =
    coefficient_indices_of_overlapping_elements(dict, gamma[findall(GridArrays.boundary_mask(gamma, domain, true))])

export nonzero_coefficients
"""
    nonzero_coefficients(dict::Dictionary1d, x::Real)

Return the coefficient indices (`UnitRange`) of the basis elements that evaluate non zero in `x`
"""
function nonzero_coefficients(dict::Dictionary1d, x::Real)
    # The init_index is the starting index of all spline elements that overlap with x
    1:length(dict)
end

# Limits of the indices of the coefficients of B that overlap with x.
# This is a tupple of (number of elements in tuple depends is equal to dimension) of two element tuples.
coefficient_index_limits_of_overlapping_elements(B::Dictionary, x::Real) =
    tuple(nonzero_coefficients(B, x))

coefficient_index_limits_of_overlapping_elements(B::TensorProductDict, x::SVector{N}) where {N} =
    [nonzero_coefficients(Bi,xi) for (Bi, xi) in zip(elements(B), x)]

# Cartesian index limits of the coefficients of B that overlap with x.
function coefficient_cartesian_index_limits_of_overlapping_elementst(B::Dictionary, x)
    range = coefficient_index_limits_of_overlapping_elements(B, x)
    CartesianIndex(range[1][1]), CartesianIndex(range[1][end])
end

function coefficient_cartesian_index_limits_of_overlapping_elementst(B::Dictionary, x::SVector{N}) where N
    ranges = coefficient_index_limits_of_overlapping_elements(B, x)
    CartesianIndex(ntuple(k->ranges[k][1] ,Val(N))), CartesianIndex(ntuple(k->ranges[k][end] ,Val(N)))
end

# Range of coefficient indices of B that overlap with the point x.
coefficient_index_range_of_overlapping_elements(B::Dictionary, x) =
    GridArrays.ModCartesianIndices(size(B), coefficient_cartesian_index_limits_of_overlapping_elementst(B, x)...)

coefficient_index_mask_of_overlapping_elements(d::Dictionary, g::AbstractGrid) =
    coefficient_index_mask_of_overlapping_elements!(BitArray(undef, size(d)), d, g)

function coefficient_index_mask_of_overlapping_elements!(mask::BitArray, B::Dictionary, g::AbstractGrid)
    fill!(mask, 0)
    for x in g, i in coefficient_index_range_of_overlapping_elements(B, x)
        mask[i] = 1
    end
    mask
end

coefficient_indices_of_overlapping_elements(dict::Dictionary, boundary::AbstractGrid) =
    findall(coefficient_index_mask_of_overlapping_elements(dict, boundary))
