module BSplineExtension
using Reexport
@reexport using CompactTranslatesDict, FrameFun, DomainSets



include("plots.jl")
include("CompactInfiniteVectors.jl")
using ..CompactInfiniteVectors

include("SparseArrayOperators.jl")
include("platforms/BSplinePlatforms.jl")
@reexport using .BSplinePlatforms
include("platforms/ndplatforms.jl")

include("BSplineExtensionSolvers.jl")
@reexport using .BSplineExtensionSolvers

FrameFun.Platforms.platform(dict::BSplineTranslatesBasis) =
    CDBSplinePlatform{coefficienttype(dict),degree(dict)}()
using CompactTranslatesDict: PeriodicEquispacedTranslates
FrameFun.Platforms.platform(dict::PeriodicEquispacedTranslates) =
    CDPETPlatform(dict)

include("AZSparse/AZSparse.jl")
@reexport using .AZSparse

include("BrainScan/BrainScan.jl")
@reexport using .BrainScan
end
