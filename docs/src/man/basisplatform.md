
# Basis platforms
## 1 dimensional B-spline platforms
Three platforms of this package generate B-spline bases as primary dictionaries, but different dual dictionaries:

```@meta
DocTestSetup = quote
    using BSplineExtension    
end
```

```jldoctest basis
julia> P1 = BSplinePlatform()
BSplinePlatform{Float64,3}()

julia> P2 = EpsBSplinePlatform()
EpsBSplinePlatform{Float64,3}()

julia> P3 = CDBSplinePlatform()
CDBSplinePlatform{Float64,3}()

julia> dictionary(P1,100)
Periodic equispaced translates of B spline of degree 3
    ↳ length = 100
    ↳ Float64 -> Float64
    ↳ support = 0.0..1.0 (Unit)
    ↳ degree = 3

julia> dictionary(P2,100)
Periodic equispaced translates of B spline of degree 3
    ↳ length = 100
    ↳ Float64 -> Float64
    ↳ support = 0.0..1.0 (Unit)
    ↳ degree = 3

julia> dictionary(P3,100)
Periodic equispaced translates of B spline of degree 3
    ↳ length = 100
    ↳ Float64 -> Float64
    ↳ support = 0.0..1.0 (Unit)
    ↳ degree = 3

julia> azdual_dict(P1,100)
Dictionary M * P

P	:	Periodic equispaced translates of B spline of degree 3
		    ↳ length = 100
		    ↳ Float64 -> Float64
		    ↳ support = 0.0..1.0 (Unit)
		    ↳ degree = 3
M	:	Multiplication by Circulant{Float64,Complex{Float64}}

julia> azdual_dict(P2,100)
Dictionary M * P

P	:	Periodic equispaced translates of B spline of degree 3
		    ↳ length = 100
		    ↳ Float64 -> Float64
		    ↳ support = 0.0..1.0 (Unit)
		    ↳ degree = 3
M	:	Multiplication by BasisFunctions.VerticalBandedMatrix{Float64}

julia> azdual_dict(P3,100)
Equispaced translates of a discrete kernel dual to B-spline
    ↳ length = 100
    ↳ Float64 -> Float64
    ↳ support = 0.0..1.0 (Unit)
    ↳ degree = 3
    ↳ m = 2
```


## 1 dimensional generic PET platforms
Two platforms of the above are generalized to `PeriodicEquispacedTranslates`:

```jldoctest basis
julia> PET = BSplineTranslatesBasis(10,3,-1,1)
Periodic equispaced translates of a periodic kernel function
    ↳ length = 10
    ↳ Float64 -> Float64
    ↳ support = -1.0..1.0

julia> P4 = PETPlatform(PET)
PETPlatform{Float64,Float64,GenericPeriodicEquispacedTranslates{Float64,Float64}}(Periodic equispaced translates of a periodic kernel function
    ↳ length = 10
    ↳ Float64 -> Float64
    ↳ support = -1.0..1.0

)

julia> P5 = CDPETPlatform(PET)
CDPETPlatform{Float64,Float64,GenericPeriodicEquispacedTranslates{Float64,Float64}}(Periodic equispaced translates of a periodic kernel function
    ↳ length = 10
    ↳ Float64 -> Float64
    ↳ support = -1.0..1.0

)

julia> dictionary(P4,100)
Periodic equispaced translates of a periodic kernel function
    ↳ length = 100
    ↳ Float64 -> Float64
    ↳ support = -1.0..1.0

julia> dictionary(P5,100)
Periodic equispaced translates of a periodic kernel function
    ↳ length = 100
    ↳ Float64 -> Float64
    ↳ support = -1.0..1.0

julia> azdual_dict(P4,100)
Dictionary M * P

P	:	Periodic equispaced translates of a periodic kernel function
		    ↳ length = 100
		    ↳ Float64 -> Float64
		    ↳ support = -1.0..1.0
M	:	Multiplication by Circulant{Float64,Complex{Float64}}

julia> azdual_dict(P5,100)
Equispaced translates of a discrete kernel dual
    ↳ Periodic equispaced translates of a periodic kernel function
      ↳ length = 100
      ↳ Float64 -> Float64
      ↳ support = -1.0..1.0
    ↳ m = 2
```

## N dimensional
```jldoctest basis
julia> P1 = NdBSplinePlatform((3,2))
ProductPlatform{2}((BSplinePlatform{Float64,3}(), BSplinePlatform{Float64,2}()))

julia> P2 = NdEpsBSplinePlatform((1,3))
ProductPlatform{2}((EpsBSplinePlatform{Float64,1}(), EpsBSplinePlatform{Float64,3}()))

julia> P3 = NdCDBSplinePlatform((3,3))
ProductPlatform{2}((CDBSplinePlatform{Float64,3}(), CDBSplinePlatform{Float64,3}()))

julia> dictionary(P1,100)
Dictionary P₂ ⊗ P₁

P₂	:	Periodic equispaced translates of B spline of degree 3
		    ↳ length = 100
		    ↳ Float64 -> Float64
		    ↳ support = 0.0..1.0 (Unit)
		    ↳ degree = 3
P₁	:	Periodic equispaced translates of B spline of degree 2
		    ↳ length = 100
		    ↳ Float64 -> Float64
		    ↳ support = 0.0..1.0 (Unit)
		    ↳ degree = 2

julia> dictionary(P2,100)
Dictionary P₁ ⊗ P₂

P₂	:	Periodic equispaced translates of B spline of degree 3
		    ↳ length = 100
		    ↳ Float64 -> Float64
		    ↳ support = 0.0..1.0 (Unit)
		    ↳ degree = 3
P₁	:	Periodic equispaced translates of B spline of degree 1
		    ↳ length = 100
		    ↳ Float64 -> Float64
		    ↳ support = 0.0..1.0 (Unit)
		    ↳ degree = 1

julia> dictionary(P3,100)
Dictionary P ⊗ P

P	:	Periodic equispaced translates of B spline of degree 3
		    ↳ length = 100
		    ↳ Float64 -> Float64
		    ↳ support = 0.0..1.0 (Unit)
		    ↳ degree = 3

julia> azdual_dict(P1,100)
Dictionary (M₁ * P₂) ⊗ (M₂ * P₁)

P₂	:	Periodic equispaced translates of B spline of degree 3
		    ↳ length = 100
		    ↳ Float64 -> Float64
		    ↳ support = 0.0..1.0 (Unit)
		    ↳ degree = 3
P₁	:	Periodic equispaced translates of B spline of degree 2
		    ↳ length = 100
		    ↳ Float64 -> Float64
		    ↳ support = 0.0..1.0 (Unit)
		    ↳ degree = 2
M₂	:	Multiplication by Circulant{Float64,Complex{Float64}}
M₁	:	Multiplication by Circulant{Float64,Complex{Float64}}

julia> azdual_dict(P2,100)
Dictionary (M₁ * P₁) ⊗ (M₂ * P₂)

P₂	:	Periodic equispaced translates of B spline of degree 3
		    ↳ length = 100
		    ↳ Float64 -> Float64
		    ↳ support = 0.0..1.0 (Unit)
		    ↳ degree = 3
P₁	:	Periodic equispaced translates of B spline of degree 1
		    ↳ length = 100
		    ↳ Float64 -> Float64
		    ↳ support = 0.0..1.0 (Unit)
		    ↳ degree = 1
M₂	:	Multiplication by BasisFunctions.VerticalBandedMatrix{Float64}
M₁	:	Multiplication by BasisFunctions.VerticalBandedMatrix{Float64}

julia> azdual_dict(P3,100;oversamplingfactor=4)
Dictionary E ⊗ E

E	:	Equispaced translates of a discrete kernel dual to B-spline
		    ↳ length = 100
		    ↳ Float64 -> Float64
		    ↳ support = 0.0..1.0 (Unit)
		    ↳ degree = 3
		    ↳ m = 2
```

# Function approximation

These duals are used to approximate functions if the `SolverStyle` option is equal to `DualStyle` (default):


[\[.pdf\]](figs/basis.pdf), [\[generated .tex\]](figs/basis.tex), [\[generated .tikz\]](figs/basis.tikz)

![](figs/basis.svg)

## Convergence

The convergence curves show no difference between the different approaches:

[\[.pdf\]](figs/convergence_basis.pdf), [\[generated .tex\]](figs/convergence_basis.tex), [\[generated .tikz\]](figs/convergence_basis.tikz)

![](figs/convergence_basis.svg)

# B-Spline BasisPlatforms Reference

```@autodocs
Modules = [BSplineExtension.BSplinePlatforms]
Pages = ["platforms/BSplinePlatforms.jl"]
```

The dual dictionary used in [`CDBSplinePlatform`](@ref) is

```@autodocs
Modules = [BSplineExtension.BSplinePlatforms]
Pages = ["platforms/DiscreteBSplineDicts.jl"]
```

and the one used in [`CDPETPlatform`](@ref) is
```@autodocs
Modules = [BSplineExtension.BSplinePlatforms.CompactPeriodicEquispacedTranslatesDuals]
Pages = ["platforms/CompactPeriodicEquispacedTranslatesDuals.jl"]
```


```@autodocs
Modules = [BSplineExtension]
Pages = ["platforms/ndplatforms.jl"]
```
