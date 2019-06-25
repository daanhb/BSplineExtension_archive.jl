
import FrameFun.BasisFunctions: string, strings, name, unsafe_eval_element, grid_evaluation_operator, size, length
using FrameFun.BasisFunctions: VerticalBandedMatrix

using InfiniteVectors: Integers, downsample, subvector
using CardinalBSplines: compact_dualperiodicbsplinesignal
string(::Integers) = "The integers (ℤ)"
using CompactTranslatesDict: PeriodicBSplineBasis

export DiscreteBSplineDict
"""
    struct DiscreteBSplineDict{T,K} <: PeriodicBSplineBasis{T,K}

Compact, discrete, periodic dual --- with respect to a `DiracCombMeasure` ---
to the B-spline basis consisting of equispaced
periodic translates of a B-spline (see `BSplineTranslatesBasis` of `CompactTranslatesDict`).

Since this basis is discrete, it can only be evaluated in an `PeriodicEquispacedGrid`.

Used in [`CDBSplinePlatform`](@ref)

# Example
```jldocs
julia> B = DiscreteBSplineDict(3,2,10)
Equispaced translates of a discrete kernel dual to B-spline
    ↳ length = 10
    ↳ Float64 -> Float64
    ↳ support = 0.0..1.0 (Unit)
    ↳ degree = 3
    ↳ m = 2



julia> A = evaluation_operator(B, PeriodicEquispacedGrid(10, support(B)))
Multiplication by BasisFunctions.VerticalBandedMatrix{Float64}



julia> Matrix(A)
10×10 Array{Float64,2}:
  0.38283     -0.17734      0.11106      0.00547352   0.0          0.0          0.0          0.00547352   0.11106     -0.17734
 -0.17734      0.38283     -0.17734      0.11106      0.00547352   0.0          0.0          0.0          0.00547352   0.11106
  0.11106     -0.17734      0.38283     -0.17734      0.11106      0.00547352   0.0          0.0          0.0          0.00547352
  0.00547352   0.11106     -0.17734      0.38283     -0.17734      0.11106      0.00547352   0.0          0.0          0.0
  0.0          0.00547352   0.11106     -0.17734      0.38283     -0.17734      0.11106      0.00547352   0.0          0.0
  0.0          0.0          0.00547352   0.11106     -0.17734      0.38283     -0.17734      0.11106      0.00547352   0.0
  0.0          0.0          0.0          0.00547352   0.11106     -0.17734      0.38283     -0.17734      0.11106      0.00547352
  0.00547352   0.0          0.0          0.0          0.00547352   0.11106     -0.17734      0.38283     -0.17734      0.11106
  0.11106      0.00547352   0.0          0.0          0.0          0.00547352   0.11106     -0.17734      0.38283     -0.17734
 -0.17734      0.11106      0.00547352   0.0          0.0          0.0          0.00547352   0.11106     -0.17734      0.38283

```
"""
struct DiscreteBSplineDict{T,K} <: PeriodicBSplineBasis{T,K}
    m      ::   Int
    n      ::   Int
    function DiscreteBSplineDict{T}(degree::Int, m::Int, n::Int) where T
        new{T,degree}(m,n)
    end
    DiscreteBSplineDict(degree::Int, m::Int, n::Int) =
        DiscreteBSplineDict{Float64}(degree, m, n)
end
size(dict::DiscreteBSplineDict) = (dict.n,)
length(dict::DiscreteBSplineDict) = dict.n

name(::DiscreteBSplineDict)= "Equispaced translates of a discrete kernel dual to B-spline"

strings(d::DiscreteBSplineDict) = (string(d),
    ("length = $(length(d))",
     "$(domaintype(d)) -> $(codomaintype(d))",
     "support = $(support(d))",
     "degree = $(degree(d))",
     "m = $(d.m)"),
     )

unsafe_eval_element(dict::DiscreteBSplineDict, i, x) =
    error("`DiscreteBSplineDict` can only be evaluated in `PeriodicEquispacedGrid`")
grid_evaluation_operator(dict::DiscreteBSplineDict, gb::GridBasis, grid::AbstractGrid) =
    error("`DiscreteBSplineDict` can only be evaluated in `PeriodicEquispacedGrid`")
grid_evaluation_operator(dict::DiscreteBSplineDict, dgs::GridBasis, grid::AbstractEquispacedGrid) =
    error("`DiscreteBSplineDict` can only be evaluated in `PeriodicEquispacedGrid`")
eval_kernel(dict::DiscreteBSplineDict, x) =
    error("`DiscreteBSplineDict` can only be evaluated in `PeriodicEquispacedGrid`") 

function grid_evaluation_operator(dict::DiscreteBSplineDict{T}, gb::GridBasis, grid::PeriodicEquispacedGrid; options...) where {T}
    @assert support(dict) ≈ support(grid)
    m_dict = dict.m
    m_frac = m_dict*length(dict) / length(grid)

    @assert round(m_frac) ≈ m_frac
    m_frac = round(Int, m_frac)


    v = downsample(compact_dualperiodicbsplinesignal(degree(dict), m_dict, m_dict*length(dict), T), m_frac)

    # Convert PeriodicInfiniteVector to VerticalBandedMatrix
    mask = 0 .==subvector(v)
    bw = findfirst(mask)-1, findlast(mask)+1
    a = [subvector(v)[bw[2]:end];subvector(v)[1:bw[1]]]
    if true # BSplineTranslatesBasis is scaled
        a ./= sqrt(length(dict))
    end


    (length(grid)<length(dict)) && (error("Not implemented."))
    shift = length(grid) / length(dict)
    (round(shift) ≈ shift) || error("not implemented")
    shift = round(Int, shift)

    # v is calculated for centered bspline. we need a shift to align with the elements of `BSplineTranslatesBasis`: m_dict*degree...



    ArrayOperator(VerticalBandedMatrix(length(grid), length(dict), a, shift, -1+bw[2]-length(grid)), dict, GridBasis{coefficienttype(dict)}(grid))
end
