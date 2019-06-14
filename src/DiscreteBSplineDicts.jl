
import FrameFun.BasisFunctions: string, strings, name, unsafe_eval_element, grid_evaluation_operator, size, length
using FrameFun.BasisFunctions: VerticalBandedMatrix

using InfiniteVectors: Integers, downsample, subvector
using CardinalBSplines: compact_dualperiodicbsplinesignal
string(::Integers) = "The integers (ℤ)"
using CompactTranslatesDict: PeriodicBSplineBasis

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
