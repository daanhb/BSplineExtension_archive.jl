module BSplineExtension
using Reexport, LinearAlgebra, StaticArrays
@reexport using CompactTranslatesDict, FrameFun, DomainSets

using FrameFun.FrameFunInterface: directsolver
import FrameFun: dictionary
import FrameFun.Platforms: dictionary, SolverStyle, measure, SamplingStyle, dualdictionary
import FrameFun.FrameFunInterface: correct_sampling_parameter

include("plots.jl")

include("platforms/AbstractBSplinePlatforms.jl")

function nonzero_coefficients(dict::CompactTranslatesDict.DiffPeriodicBSplineBasis, x::Real)
    w = (degree(dict)+1)
    os = iseven(degree(dict)) ? -.5 : 0
    tuple(floor(Int, 1 + os + x/step(dict)) .+ (-(w>>1)+1:((w+1)>>1)))
end

include("BSplineExtensionSolver.jl")

include("platforms/ndplatforms.jl")


end
