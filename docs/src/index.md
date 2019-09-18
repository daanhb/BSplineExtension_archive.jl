
# BSplineExtension.jl documentation

*A Julia package for B-spline extension approximations*


## Installation

`BSplineExtension.jl` is not added to the Julia General registry, it is however in the FrameFun registry.

### Recomanded
For Julia 1.1 or higher, you can add the FrameFun registry to the list of your registries. Then add `BSplineExtension.jl`.
From the Julia REPL, type `]` to enter Pkg mode and run

```julia
pkg> registry add https://github.com/FrameFunVC/FrameFunRegistry
pkg> add BSplineExtension
```

### Legacy
In Julia 1.0, the package and its dependencies can be installed by cloning their git repository. From the Julia REPL, type `]` to enter Pkg mode and run

```julia
pkg> add https://github.com/FrameFunVC/InfiniteVectors.jl
pkg> add https://github.com/FrameFunVC/CardinalBSplines.jl
pkg> add https://github.com/JuliaApproximation/BasisFunctions.jl
pkg> add https://github.com/FrameFunVC/CompactTranslatesDict.jl
pkg> add https://github.com/JuliaApproximation/FrameFun.jl
```

## Manual
This module consists of multiple platforms. Three of them are basis platforms. Each of them provide a `BSplineTranslatesBasis` if asked for a dictionary; but they all provide a different kind of dual dictionary, see [Basis platforms](@ref). The module also extends the `FrameFun` framework to 1 dimensional B-spline extensions, see [1D B-spline extension](@ref).
