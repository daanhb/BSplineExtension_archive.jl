module BrainPlatforms

using FrameFun.Platforms, GridArrays, ...AZSparse, FrameFun.ExtensionFrames, FrameFun.BasisFunctions,
    DomainSets, ..BrainScanGrids
using DomainSets: WrappedDomain
import FrameFun.Platforms: platform, SolverStyle, SamplingStyle, measure, dictionary, dualdictionary
import FrameFun.FrameFunInterface: platform_grid
import FrameFun.ExtensionFrames: support


export BrainPlatform
struct BrainPlatform <: FramePlatform
    basisplatform::Platform
    grid::AbstractGrid
    support::Domain
    opts::Base.Iterators.Pairs
    function BrainPlatform(basisplatform::Platform;supportgrid=nothing,domain=nothing,opts...)
        if supportgrid!=nothing || domain!=nothing
            @warn "brain grid uses a default domain and supportgrid"
        end
        grid = braingrid(;domain=support(dictionary(basisplatform,(1,1,1))), opts...)
        new(basisplatform, grid, WrappedDomain(grid), opts)
    end
end

support(platform::BrainPlatform) =
    error("better not to use the support of a BrainPlatform")
platform_grid(ss::GridStyle, platform::BrainPlatform, param; options...) =
    platform.grid



SolverStyle(dict::BrainPlatform, ::SamplingStyle) = AZSparseStyle()

SamplingStyle(p::BrainPlatform) = GridStyle()

dictionary(p::BrainPlatform, n) =
    extensionframe(p.support, dictionary(p.basisplatform, n))

measure(platform::BrainPlatform) =
    discretemeasure(platform.grid)

dualdictionary(platform::BrainPlatform, param, measure::Measure; options...) =
   extensionframe(dualdictionary(platform.basisplatform, param, supermeasure(measure); options...), platform.support)




end
