
# Basis platforms
## 1 dimensional
Three platforms of this package generate B-spline bases as primary dictionaries, but different dual dictionaries:

```@meta
DocTestSetup = quote
    using BSplineExtension    
end
```

```@setup basis
using PGFPlotsX, BSplineExtension
using BSplineExtension.FrameFun: abserror
P1 = BSplinePlatform()
P2 = EpsBSplinePlatform()
P3 = CDBSplinePlatform()
savefigs = (figname, obj) -> begin
    pgfsave(figname * ".pdf", obj)
    pgfsave(figname * ".tex", obj);
    pgfsave(figname * ".tikz", obj;include_preamble=false);
    run(`pdf2svg $(figname * ".pdf") $(figname * ".svg")`)
    return nothing
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

## N dimensional
```jldoctest basis
julia> P1 = NdBSplinePlatform((3,2))
ProductPlatform{2}((BSplinePlatform{Float64,3}(), BSplinePlatform{Float64,2}()))

julia> P2 = NdEpsBSplinePlatform((1,3))
ProductPlatform{2}((EpsBSplinePlatform{Float64,1}(), EpsBSplinePlatform{Float64,3}()))

julia> P3 = NdCDBSplinePlatform((3,3))
ProductPlatform{2}((CDBSplinePlatform{Float64,3}(), CDBSplinePlatform{Float64,3}()))

julia> dictionary(P1,100)
Dictionary P₁ ⊗ P₂

P₂	:	Periodic equispaced translates of B spline of degree 2
		    ↳ length = 100
		    ↳ Float64 -> Float64
		    ↳ support = 0.0..1.0 (Unit)
		    ↳ degree = 2
P₁	:	Periodic equispaced translates of B spline of degree 3
		    ↳ length = 100
		    ↳ Float64 -> Float64
		    ↳ support = 0.0..1.0 (Unit)
		    ↳ degree = 3

julia> dictionary(P2,100)
Dictionary P₂ ⊗ P₁

P₂	:	Periodic equispaced translates of B spline of degree 1
		    ↳ length = 100
		    ↳ Float64 -> Float64
		    ↳ support = 0.0..1.0 (Unit)
		    ↳ degree = 1
P₁	:	Periodic equispaced translates of B spline of degree 3
		    ↳ length = 100
		    ↳ Float64 -> Float64
		    ↳ support = 0.0..1.0 (Unit)
		    ↳ degree = 3

julia> dictionary(P3,100)
Dictionary P ⊗ P

P	:	Periodic equispaced translates of B spline of degree 3
		    ↳ length = 100
		    ↳ Float64 -> Float64
		    ↳ support = 0.0..1.0 (Unit)
		    ↳ degree = 3

julia> azdual_dict(P1,100)
Dictionary (M₂ * P₁) ⊗ (M₁ * P₂)

P₂	:	Periodic equispaced translates of B spline of degree 2
		    ↳ length = 100
		    ↳ Float64 -> Float64
		    ↳ support = 0.0..1.0 (Unit)
		    ↳ degree = 2
P₁	:	Periodic equispaced translates of B spline of degree 3
		    ↳ length = 100
		    ↳ Float64 -> Float64
		    ↳ support = 0.0..1.0 (Unit)
		    ↳ degree = 3
M₂	:	Multiplication by Circulant{Float64,Complex{Float64}}
M₁	:	Multiplication by Circulant{Float64,Complex{Float64}}

julia> azdual_dict(P2,100)
Dictionary (M₁ * P₂) ⊗ (M₂ * P₁)

P₂	:	Periodic equispaced translates of B spline of degree 1
		    ↳ length = 100
		    ↳ Float64 -> Float64
		    ↳ support = 0.0..1.0 (Unit)
		    ↳ degree = 1
P₁	:	Periodic equispaced translates of B spline of degree 3
		    ↳ length = 100
		    ↳ Float64 -> Float64
		    ↳ support = 0.0..1.0 (Unit)
		    ↳ degree = 3
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

```@example basis
f = x->exp(cos(10pi*x));N=60;
F1 = Fun(f, P1, N)
F2 = Fun(f, P2, N)
F3 = Fun(f, P3, N)
@pgf GroupPlot({group_style = {group_size="2 by 1",},},
        {},
        PlotInc(F1;n=600),
        LegendEntry("P1"),
        PlotInc(F2;n=600),
        LegendEntry("P2"),
        PlotInc(F3;n=600),
        LegendEntry("P3"),
        {ymode="log"},
        PlotInc(f, F1;n=600),
        PlotInc(f, F2;n=600),
        PlotInc(f, F3;n=600));
savefigs("basis", ans) # hide
```

[\[.pdf\]](basis.pdf), [\[generated .tex\]](basis.tex), [\[generated .tikz\]](basis.tikz)

![](basis.svg)

## Convergence

The convergence curves show no difference between the different approaches:
```@example basis
Ns = 10*[1<<k for k in 0:11] # hide
err = zeros(3, length(Ns)) # hide
for (i,P) in enumerate((P1,P2,P3))
    for (j,N) in enumerate(Ns)
        F = Fun(f, P, N)
        err[i,j] = abserror(f, F)
    end
end
@pgf Axis({ymode="log",xmode="log"},
        PlotInc(Table(Ns,err[1,:])),
        LegendEntry("P1"),
        PlotInc(Table(Ns,err[2,:])),
        LegendEntry("P2"),
        PlotInc(Table(Ns,err[3,:])),
        LegendEntry("P3"),
        )
savefigs("convergence_basis", ans) # hide
```

[\[.pdf\]](convergence_basis.pdf), [\[generated .tex\]](convergence_basis.tex), [\[generated .tikz\]](convergence_basis.tikz)

![](convergence_basis.svg)

# B-Spline BasisPlatforms Reference

```@autodocs
Modules = [BSplineExtension]
Pages = ["platforms/AbstractBSplinePlatforms.jl"]
```

The dual dictionary used in [`CDBSplinePlatform`](@ref) is

```@autodocs
Modules = [BSplineExtension]
Pages = ["platforms/DiscreteBSplineDicts.jl"]
```

```@autodocs
Modules = [BSplineExtension]
Pages = ["platforms/ndplatforms.jl"]
```
