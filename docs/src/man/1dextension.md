
# 1D B-spline extension
```@setup 1dframe
using PGFPlotsX, BSplineExtension, DomainSets, LaTeXStrings
P1 = ExtensionFramePlatform(BSplinePlatform(),0.0..0.5)
P2 = ExtensionFramePlatform(EpsBSplinePlatform(),0.0..0.5)
P3 = ExtensionFramePlatform(CDBSplinePlatform(),0.0..0.5)
savefigs = (figname, obj) -> begin
    pgfsave(figname * ".pdf", obj)
    run(`pdf2svg $(figname * ".pdf") $(figname * ".svg")`)
    pgfsave(figname * ".tex", obj);
    return nothing
end
```

```@meta
DocTestSetup = quote
    using BSplineExtension, DomainSets  
end
```

The package provides a fast solver for a B-spline extensions such as, e.g.,
```jldoctest 1dframe
julia> P1 = ExtensionFramePlatform(BSplinePlatform(),0.0..0.5)
ExtensionFramePlatform(BSplinePlatform{Float64,3}(), 0.0..0.5)
julia> P2 = ExtensionFramePlatform(EpsBSplinePlatform(),0.0..0.5)
ExtensionFramePlatform(EpsBSplinePlatform{Float64,3}(), 0.0..0.5)
julia> P3 = ExtensionFramePlatform(CDBSplinePlatform(),0.0..0.5)
ExtensionFramePlatform(CDBSplinePlatform{Float64,3}(), 0.0..0.5)
```

## The AZ algorithm
The AZ algorithm (reference needed) provides for a fast approximation algorithm
even though the system that needs to be solved is rectangular, and highly ill conditioned

The AZ algorithm consists of three steps for solving a system `Ax=b`
- Solve `(A-AZ'A)x1 = (I-AZ')b`
- `x2 = Z'(b-Ax1)`
- `x = x1 + x2`

Z can be any matrix, but the algorithm is most efficient if `A` and `Z'` can be applied fast, and
`A-AZ'A` is of low rank.

For the platforms above these operators are readily available using `AZ_A` and `AZ_Zt`:
```jldoctest 1dframe
julia> N = 200;


julia> AZ_A(P1,N)
Operator D * R[1:401 → 𝕀] * M

R	:	Restriction of coefficients to subset
M	:	Multiplication by BasisFunctions.VerticalBandedMatrix{Float64}
D	:	Diagonal operator with element type Float64
		    ↳ [1.0, 1.0, 1.0  …  1.0, 1.0, 1.0]

julia> AZ_Zt(P1,N)
Operator M₁ * M₂ * E[ 𝕀 → 1:401] * D

M₂	:	Multiplication by BasisFunctions.HorizontalBandedMatrix{Float64}
M₁	:	Multiplication by Circulant{Float64,Complex{Float64}}
E	:	Extending coefficients by zero padding
D	:	Diagonal operator with element type Float64
		    ↳ [1.0, 1.0, 1.0  …  1.0, 1.0, 1.0]

julia> AZ_Zt(P2,N)
Operator M₂ * M₁ * E[ 𝕀 → 1:401] * D

M₂	:	Multiplication by BasisFunctions.HorizontalBandedMatrix{Float64}
M₁	:	Multiplication by BasisFunctions.HorizontalBandedMatrix{Float64}
E	:	Extending coefficients by zero padding
D	:	Diagonal operator with element type Float64
		    ↳ [1.0, 1.0, 1.0  …  1.0, 1.0, 1.0]


julia> AZ_Zt(P3,N)
Operator M * E[ 𝕀 → 1:401] * D

M	:	Multiplication by BasisFunctions.HorizontalBandedMatrix{Float64}
E	:	Extending coefficients by zero padding
D	:	Diagonal operator with element type Float64
		    ↳ [1.0, 1.0, 1.0  …  1.0, 1.0, 1.0]


```

## Sparseness of `A` and `Z'`
The `FrameFun` package contains fast solvers for solving the low-rank system of
the first AZ step. However, it does not use the (ϵ-)sparse (with sparse: a lot of zero elements)
structure of `A` and `Zt`
which is a consequence of the compact nature of a B-spline basis.
For all platforms `AZ_A` returns the same operator. As is seen above `A` is a banded matrix, hence sparse.
For `P1` is `Z'` is the product of a banded matrix and a circulant matrix.
The circulant matrix has exponentially decreasing elements in the first column, hence ϵ-sparse.
Only for `P2` and `P3`, `Z'` is actually be sparse.

# BSplineExtensionSolver
The [`BSplineExtensionSolver`](@ref) of this package uses the compactness of the B-spline
basis and by consequence the sparseness of `A` and `Z'`.
First, the compactness of the B-spline basis introduces a lot of zero columns in `A-AZ'A`.
Therefore, the BSplineExtensionSolver removes the zero columns. It knows of the location of
the non zero columns by using [`nonzero_coefficients`](@ref).
Secondly, since also the dual dictionary is (ϵ-)compactness (see [Basis platforms](@ref) for examples of these dual dictionaries),
`A-AZ'A` will contain some (ϵ-)zero rows. These are also removed by the `BSplineExtensionSolver`, where `crop_tol` is a parameter that
indicates how large ϵ should be. The resulting system is of size [`truncated_size`](@ref).

Below we illustrate how the `truncated_size` of a `BSplineExtensionSolver` changes for varying
platform, B-spline order (`m`), `crop_tol` and degrees of freedom (`N`).

## Column truncation
The column truncation of `A-AZ'A` depends on the order of the B-spline (and the domain) only.
This is seen in the figure below, which shows `truncated_size(S)[2]` with `S=BSplineExtensionSolver(A-A*Zt*A)`.
```@example 1dframe
Ns = 20:20:300 # hide
ds = 1:4 # hide
PLATFORMs = (BSplinePlatform, EpsBSplinePlatform, CDBSplinePlatform) # hide
crop_tols = 10.0.^(-16.:6.:-10.) # hide
colsizes = Array{Int}(undef, length(PLATFORMs), length(ds), length(Ns), length(crop_tols)) # hide
rowsizes = similar(colsizes) # hide
for (i,PLATFORM) in enumerate(PLATFORMs), (j,d) in enumerate(ds), (k,N) in enumerate(Ns), (l,crop_tol) in enumerate(crop_tols) # hide
    P = ExtensionFramePlatform(PLATFORM(d), 0.0..0.5); # hide
    plunge = plungeoperator(P,N;L=4N); A = AZ_A(P,N;L=4N); Zt = AZ_Zt(P,N;L=4N); # hide
    M = plunge*A; # hide
    colsizes[i,j,k,l], rowsizes[i,j,k,l]  = truncated_size(BSplineExtensionSolver(M; crop_tol=crop_tol)) # hide
end # hide
A = [];for (i,PLATFORM) in enumerate(PLATFORMs) # hide
    push!(A,@pgf {legend_pos="north west",xlabel="N",title=["BSplinePlatform","EpsBSplinePlatform","CDBSplinePlatform"][i]}) # hide
    for (j,d) in enumerate(ds) # hide
        push!(A, @pgf PlotInc({color=["blue", "red", "brown", "black"][j],mark=["o", "square", "diamond", "x"][j],mark_options="solid"}, Table(Ns, rowsizes[i,j,:,1]))) # hide
        i==1 && push!(A, @pgf LegendEntry("m=$d")) # hide
    end # hide
end # hide
@pgf PGFPlotsX.GroupPlot({ymin=0,group_style={group_size="3 by 1",},}, # hide
    A...) # hide
savefigs("truncated_size_1", ans) # hide
```

[\[.pdf\]](truncated_size_1.pdf), [\[generated .tex\]](truncated_size_1.tex)

![](truncated_size_1.svg)

## Row truncation
The row truncation of `A-AZ'A` depends on all of the parameters.
This is seen in the figure below, which shows `truncated_size(S)[1]` with
`S=BSplineExtensionSolver(A-A*Zt*A;crop_tol=crop_tol)`
for `crop_tol` in (1e-16, 1e-10).

The truncated size reaches a limit for `N` growing large. This limit is much smaller using `CDBSplinePlatform`.
The `EpsBSplinePlatform` is a bit more robust  than the `BSplinePlatform` in generating a low truncated size
since it enforces the dual dictionary to be
compact, while the dual dictionaries of the latter are only ϵ-compact.

```@example 1dframe
A = [];for (i,PLATFORM) in enumerate(PLATFORMs) # hide
    push!(A,@pgf {xlabel="N",legend_pos="north west",title=["BSplinePlatform","EpsBSplinePlatform","CDBSplinePlatform"][i]}) # hide
    for (j,d) in enumerate(ds) # hide
        opts = @pgf {color=["blue", "red", "brown", "black"][j],mark=["o", "square", "diamond", "x"][j],mark_options="solid"} # hide
        for (l,crop_tol) in enumerate(crop_tols) # hide
            push!(A, @pgf Plot({opts..., (@pgf {solid}, {dashed})[l]...}, Table(Ns, colsizes[i,j,:,l]))) # hide
            i==1 && push!(A, @pgf LegendEntry("m=$d($(crop_tol))")) # hide
        end # hide
    end # hide
end # hide
@pgf PGFPlotsX.GroupPlot({ymin=0,group_style={group_size="3 by 1",},}, # hide
    A...) # hide
savefigs("truncated_size_2", ans) # hide
```

[\[.pdf\]](truncated_size_2.pdf), [\[generated .tex\]](truncated_size_2.tex)

![](truncated_size_2.svg)

# B-spline extension approximation
The [`BSplineExtensionSolver`](@ref) is just a means to an end. It is used in the
first step of the AZ algorithm. In this section we use the solver to approximate
the exponential function on the interval [0,0.5] using a B-spline basis of order
`m` on the interval [0,1]. First we show the convergence results, then the time complexity
of the approximation algorithm.

## Errors
In the figure below, which shows the uniform error of approximating a 1 dimensional analytic function (details are in the introduction of this section),
 we see that convergence is algebraic ($$\mathcal O(N^{-m})$$).
This is common behaviour for approximating with splines.

```@example 1dframe
PLATFORMs = (EpsBSplinePlatform, BSplinePlatform, CDBSplinePlatform) # hide
Ns1 = [1<<k for k in 4:10] # hide
Ns2 = [1<<k for k in 9:16] # hide
ds = 1:4 # hide
errors = Array{Float64}(undef, length(PLATFORMs), length(ds), length(Ns1)) # hide
timings = Array{Float64}(undef, length(PLATFORMs), length(ds), length(Ns2)) # hide
for (i,PLATFORM) in enumerate(PLATFORMs), (j,d) in enumerate(ds), (k,N) in enumerate(Ns1) #hide
    P = ExtensionFramePlatform(PLATFORM(d), 0.0..0.5) #hide
    F,_ = @timed Fun(exp, P, N;L=4N, REG=BSplineExtension.BSplineExtensionSolver, crop=true, crop_tol=1e-10) #hide
    errors[i,j,k] = abserror(exp, F) #hide
end #hide
for (i,PLATFORM) in enumerate(PLATFORMs), (j,d) in enumerate(ds), (k,N) in enumerate(Ns2) #hide
    P = ExtensionFramePlatform(PLATFORM(d), 0.0..0.5); #hide
    _,timings[i,j,k],_ = @timed Fun(exp, P, N;L=4N, REG=BSplineExtension.BSplineExtensionSolver,  crop=true, crop_tol=1e-10) #hide
end # hide
A = [];for (i,PLATFORM) in enumerate(PLATFORMs) # hide
    push!(A,@pgf {xmode="log",ymode="log",xlabel="N",legend_pos="north west",title=["BSplinePlatform","EpsBSplinePlatform","CDBSplinePlatform"][i]}) # hide
    for (j,d) in enumerate(ds) # hide
        opts = @pgf {color=["blue", "red", "brown", "black"][j],mark=["o", "square", "diamond", "x"][j],mark_options="solid"} # hide
        push!(A, @pgf Plot(opts, Table(Ns1, errors[i,j,:]))) # hide
        i==1 && push!(A, @pgf LegendEntry("m=$d")) # hide
    end # hide
end # hide
@pgf PGFPlotsX.GroupPlot({ymin=0,group_style={group_size="3 by 1",},}, # hide
    A...) # hide
savefigs("1derrors", ans) # hide
```

[\[.pdf\]](1derrors.pdf), [\[generated .tex\]](1derrors.tex)

![](1derrors.svg)
## Timings
For all platforms, the most costly part is the evaluating operator-vector multiply ``(A-AZ'A)x` since
solving the system is constant in `N` (it is $$\mathcal O(d^3)$$).
For the first platform the operator-vector multiply time complexity is $$\mathcal O(N\log N)$$. For the
last two the multiply time complexity is $$\mathcal O(N)$$. However, the construction complexity of `Z'` for the first two
platforms is also $$\mathcal O(N\log N)$$.  

The quasi-linear behaviour of the method is confirmed in the figure below that
show the timings for approximating a 1 dimensional analytic function
    (details are in the introduction of this section).

```@example 1dframe
A = [];for (i,PLATFORM) in enumerate(PLATFORMs) # hide
    push!(A,@pgf {xmode="log",ymode="log",xlabel="N",legend_pos="north west",title=["BSplinePlatform","EpsBSplinePlatform","CDBSplinePlatform"][i]}) # hide
    for (j,d) in enumerate(ds) # hide
        opts = @pgf {color=["blue", "red", "brown", "black"][j],mark=["o", "square", "diamond", "x"][j],mark_options="solid"} # hide
        push!(A, @pgf Plot(opts, Table(Ns2, timings[i,j,:]))) # hide
        i==1 && push!(A, @pgf LegendEntry("m=$d")) # hide
    end # hide
    push!(A, @pgf Plot({color="black",dashed},Table(Ns2,5e-6Ns2))) # hide
    i==1 && push!(A, LegendEntry(L"\mathcal O(N)")) # hide
end # hide
@pgf PGFPlotsX.GroupPlot({ymin=0,group_style={group_size="3 by 1",},}, # hide
    A...) # hide
savefigs("1dtimings", ans) # hide
```

[\[.pdf\]](1dtimings.pdf), [\[generated .tex\]](1dtimings.tex)

![](1dtimings.svg)



# B-spline extension Reference

```@autodocs
Modules = [BSplineExtension]
Pages = ["nonzero_cols.jl"]
```

```@autodocs
Modules = [BSplineExtension]
Pages = ["BSplineExtensionSolver.jl"]
```
