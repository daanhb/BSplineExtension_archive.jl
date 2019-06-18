module BSplineExtension
using Reexport, LinearAlgebra, StaticArrays
@reexport using CompactTranslatesDict, FrameFun, DomainSets

using FrameFun: BasisPlatform, FramePlatform, ExtensionFramePlatform, Measure
import FrameFun: dictionary, first_parameters, SolverStyle, measure, azdual_dict, SamplingStyle,
    SolverStyle, oversampling_grid, deduce_samplingparameter

export ExtensionFramePlatform

include("plots.jl")

include("basisplatforms.jl")

function nonzero_coefficients(dict::CompactTranslatesDict.DiffPeriodicBSplineBasis, x::Real)
    w = (degree(dict)+1)
    os = iseven(degree(dict)) ? -.5 : 0
    floor(Int, 1 + os + x/step(dict)) .+ (-(w>>1)+1:((w+1)>>1))
end

include("BSplineExtensionSolver.jl")


end
