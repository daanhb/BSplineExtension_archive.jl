
export SparseCDBSplinePlatform
"""
    struct SparseCDBSplinePlatform{T,D} <: AbstractBSplinePlatform{T,D}

A platform of equispaced periodic translates of B-spline of a given order.
Their duals dictionaries are compact, but discrete, i.e.,
their values are known in an `PeriodicEquispacedGrid`

See also: [`DiscreteBSplineDict`](@ref), [`BSplinePlatform`](@ref), [`EpsBSplinePlatform`](@ref)

# Example
```jldocs
julia> P = CDBSplinePlatform(5)
CDBSplinePlatform{Float64,5}()

julia> d1 = dictionary(P,20)
Periodic equispaced translates of B spline of degree 5
    ↳ length = 20
    ↳ Float64 -> Float64
    ↳ support = 0.0..1.0 (Unit)
    ↳ degree = 5



julia> d2 = azdual_dict(P,20)
Equispaced translates of a discrete kernel dual to B-spline
    ↳ length = 20
    ↳ Float64 -> Float64
    ↳ support = 0.0..1.0 (Unit)
    ↳ degree = 5
    ↳ m = 2



julia> g2 = mixedgramoperator(d1, d2, discretemeasure(sampling_grid(P,20)))
Operator M₂ * M₁

M₂  :   Multiplication by BasisFunctions.HorizontalBandedMatrix{Float64}
M₁  :   Multiplication by BasisFunctions.VerticalBandedMatrix{Float64}


julia> IdentityOperator(d1)≈g2
true
```
"""
struct SparseCDBSplinePlatform{T,D} <: AbstractBSplinePlatform{T,D}
    platform    :: CDBSplinePlatform{T,D}
end

SparseCDBSplinePlatform(degree::Int=3) = SparseCDBSplinePlatform{Float64,degree}(CDBSplinePlatform(degree))

MacroTools.@forward SparseCDBSplinePlatform.platform SamplingStyle
dualdictionary(platform::SparseCDBSplinePlatform, param, measure::UniformDiracCombMeasure;
        options...) =
    dualdictionary(platform.platform, param, measure; options...)

import FrameFun.FrameFunInterface: discretization, dualdiscretization
discretization(ss::SamplingStyle, platform::SparseCDBSplinePlatform, param, S; options...) =
    sparse(discretization(ss, platform.platform, param, S; options...))
dualdiscretization(ss::SamplingStyle, platform::SparseCDBSplinePlatform, param, S, dualdict; options...) =
    sparse(dualdiscretization(ss, platform.platform, param, S, dualdict; options...))
