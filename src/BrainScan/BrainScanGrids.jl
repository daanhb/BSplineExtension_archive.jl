module BrainScanGrids
using GridArrays, DomainSets, FrameFun.BasisFunctions
using GridArrays: MaskedGrid
using DomainSets: WrappedDomain


import FrameFun.BasisFunctions: grid
grid(A::AbstractArray{T,N};domain=UnitInterval{T}()^N, supportgrid=nothing, opts...) where{T,N} =
    supportgrid==nothing ?
        grid(A, domain; opts...) :
        grid(A, supportgrid; opts...)

grid(A::AbstractArray, domain::Domain;periodic=true, opts...) =
    periodic ?
        grid(A, ProductGrid(map(PeriodicEquispacedGrid, size(A), infimum(domain), supremum(domain))...); opts...) :
        grid(A, ProductGrid(map(EquispacedGrid, size(A), infimum(domain), supremum(domain))...); opts...)

grid(A::AbstractArray, g::AbstractGrid; isin=k->!isnan(k), opts...) =
    MaskedGrid(g, convert(BitArray, isin.(A) ), WrappedDomain(g))
end
