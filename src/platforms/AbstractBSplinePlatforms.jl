abstract type AbstractBSplinePlatform{T,D} <: BasisPlatform end
dictionary(p::AbstractBSplinePlatform{T,D}, n::Int) where {T,D} = BSplineTranslatesBasis{T,D}(n)
SolverStyle(p::AbstractBSplinePlatform, ::SamplingStyle) = DualStyle()
correct_sampling_parameter(::AbstractBSplinePlatform, param, L; options...) = error()
correct_sampling_parameter(::AbstractBSplinePlatform, param::Int, L::Int; options...) =
    (round(Int, L/param) * param)


export BSplinePlatform
"""
    struct BSplinePlatform{T,D} <: AbstractBSplinePlatform{T,D}

A platform of equispaced periodic translates of B-spline of a given order.
The dual dictionaries are determined by the inverse of the Gram matrix.


See also: [`EpsBSplinePlatform`](@ref), [`CDBSplinePlatform`](@ref)

# Example
```jldocs
julia> P = BSplinePlatform()
BSplinePlatform{Float64,3}()

julia> g = sampling_grid(P,10)
10-element PeriodicEquispacedGrid{Float64}:
 0.0
 0.09999999999999999
 0.19999999999999998
 0.3
 ⋮
 0.7000000000000001
 0.7999999999999999
 0.9

julia> d1 = dictionary(P,10)
Periodic equispaced translates of B spline of degree 3
    ↳ length = 10
    ↳ Float64 -> Float64
    ↳ support = 0.0..1.0 (Unit)



julia> d2 = azdual_dict(P,10)
Dictionary M * P

P   :   Periodic equispaced translates of B spline of degree 3
            ↳ length = 10
            ↳ Float64 -> Float64
            ↳ support = 0.0..1.0 (Unit)
M   :   Multiplication by Circulant{Float64,Complex{Float64}}


julia> mixedgramoperator(d1,d2,discretemeasure(g))≈IdentityOperator(d1,d2)
true
```
"""
struct BSplinePlatform{T,D} <: AbstractBSplinePlatform{T,D}
end

BSplinePlatform(degree::Int=3) = BSplinePlatform{Float64,degree}()


export EpsBSplinePlatform
"""
    struct EpsBSplinePlatform{T,D} <: AbstractBSplinePlatform{T,D}

A platform of equispaced periodic translates of B-spline of a given order.
The dual dictionaries are determined by the inverse of the Gram matrix. This Gram
matrix has exponentially decaying elements. This platform truncates elemens lower
`threshold` which results in a banded Gram.

See also: [`BSplinePlatform`](@ref), [`CDBSplinePlatform`](@ref)

# Example
```jldocs
julia> P = EpsBSplinePlatform()
EpsBSplinePlatform{Float64,3}()

julia> d1 = dictionary(P,1000)
Periodic equispaced translates of B spline of degree 3
    ↳ length = 1000
    ↳ Float64 -> Float64
    ↳ support = 0.0..1.0 (Unit)



julia> d2 = azdual_dict(P,1000;threshold=1e-4)
Dictionary M * P

P   :   Periodic equispaced translates of B spline of degree 3
            ↳ length = 1000
            ↳ Float64 -> Float64
            ↳ support = 0.0..1.0 (Unit)
M   :   Multiplication by BasisFunctions.VerticalBandedMatrix{Float64}


julia> g2 = mixedgramoperator(d1, d2, discretemeasure(sampling_grid(P,1000)))
Operator M₂ * M₁

M₂  :   Multiplication by Circulant{Float64,Complex{Float64}}
M₁  :   Multiplication by BasisFunctions.VerticalBandedMatrix{Float64}


julia> norm(IdentityOperator(d1)-g2)
0.0004146509264381403
```
"""
struct EpsBSplinePlatform{T,D} <: AbstractBSplinePlatform{T,D}
end

EpsBSplinePlatform(degree::Int=3) = EpsBSplinePlatform{Float64,degree}()

function dualdictionary(platform::EpsBSplinePlatform, param, measure::Measure;
        threshold=nothing, options...)
    dual = gramdual(dictionary(platform, param), measure;options...)
    op = BasisFunctions.operator(dual)
    @assert op isa CirculantOperator

    # Get the abs of the first row to determine the bandwidth
    e = zeros(src(op)); e[1] = 1
    column = op*e
    # bandwidth determined by threshold
    (threshold==nothing) && (threshold=FrameFun.default_threshold(op))
    # we want \|I-G (G̃+E)\|< threshold, therefore \|E\| < threshold / \|G\|
    mask = abs.(column) .< threshold / norm(inv(op).A.vcvr_dft,Inf)
    bw = (findfirst(mask), findlast(mask))
    if bw[1] == nothing
        return dual
    end
    a = [column[bw[2]:end];column[1:bw[1]]]
    # Create bandlimited operator (which is also circulant)
    op_replace = VerticalBandedOperator(src(op), dest(op), a, 1, -1+bw[2]-length(dual))
    op_replace*src(op)
end

export CDBSplinePlatform
"""
    struct CDBSplinePlatform{T,D} <: AbstractBSplinePlatform{T,D}

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
struct CDBSplinePlatform{T,D} <: AbstractBSplinePlatform{T,D}
end

CDBSplinePlatform(degree::Int=3) = CDBSplinePlatform{Float64,degree}()

SamplingStyle(::CDBSplinePlatform) = OversamplingStyle()
include("DiscreteBSplineDicts.jl")

dualdictionary(platform::CDBSplinePlatform, param, measure::Measure; options...) =
    error("No azdual_dict for `CDBSplinePlatform` and $(typeof(measure))")

function dualdictionary(platform::CDBSplinePlatform, param, measure::UniformDiracCombMeasure;
        options...)
    dict = dictionary(platform, param)
    g = grid(measure)
    @assert support(dict) ≈ support(g)
    m = length(g) / length(dict)
    @assert round(Int,m) ≈ m
    m = round(Int, m)
    if m == 1
        @warn "No compact dual possible, try oversampling"
        return dual(dict, measure; options...)
    else
        if g isa PeriodicEquispacedGrid
            return DiscreteBSplineDict(degree(dict), m, length(dict))
        else
            return DiscreteBSplineDict(degree(dict), m, length(dict))
        end
    end
end



nonzero_azplungematrix_cols(platform::AbstractBSplinePlatform, param; dict=dictionary(platform, param), opts...) =
    _nonzero_azplungematrix_cols(dict, sampling_grid(platform, param; options...); options...)
