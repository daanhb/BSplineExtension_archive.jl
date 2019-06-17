
# BSplineExtension.jl documentation

*A Julia package for B-spline extension approximations*


## Installation

`BSplineExtension.jl` is not added to the Julia General registry, it is however in the FrameFun registry.

### Recomanded
For Julia 1.1 or higher, you can add the FrameFun registry to the list of your registries. Then add `BSplineExtension.jl`.
From the Julia REPL, type `]` to enter Pkg mode and run

```julia
pkg> add https://github.com/vincentcp/FrameFunRegistry
pkg> add BSplineExtension
```

### Legacy
In Julia 1.0, the package and its dependencies can be installed by cloning their git repository. From the Julia REPL, type `]` to enter Pkg mode and run

```julia
pkg> add https://github.com/vincentcp/InfiniteVectors.jl
pkg> add https://github.com/vincentcp/CardinalBSplines.jl
pkg> add https://github.com/vincentcp/CompactTranslatesDict.jl
pkg> add https://github.com/JuliaApproximation/FrameFun.jl
```

## Manual
This module consists of multiple platforms.
