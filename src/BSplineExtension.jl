module BSplineExtension
using Reexport
@reexport using CompactTranslatesDict, FrameFun, DomainSets



include("plots.jl")
include("SparseArrayOperators.jl")
include("platforms/BSplinePlatforms.jl")
@reexport using .BSplinePlatforms
include("platforms/ndplatforms.jl")

include("BSplineExtensionSolvers.jl")
@reexport using .BSplineExtensionSolvers

function BSplineExtension.nonzero_coefficients(dict::CompactTranslatesDict.DiffPeriodicBSplineBasis, x::Real)
    w = (degree(dict)+1)
    os = iseven(degree(dict)) ? -.5 : 0
    tuple(floor(Int, 1 + os + x/step(dict)) .+ (-(w>>1)+1:((w+1)>>1)))
end

FrameFun.Platforms.platform(dict::BSplineTranslatesBasis) =
    CDBSplinePlatform{coefficienttype(dict),degree(dict)}()
using CompactTranslatesDict: PeriodicEquispacedTranslates
FrameFun.Platforms.platform(dict::PeriodicEquispacedTranslates) =
    CDPETPlatform(dict)

include("azsparse.jl")
end
