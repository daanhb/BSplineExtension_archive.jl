# ND B-spline extension
```@setup ndframe
using PGFPlotsX, BSplineExtension, DomainSets, LaTeXStrings
P1 = ExtensionFramePlatform(NdBSplinePlatform((3,3)),(0.0..0.5)^2)
P2 = ExtensionFramePlatform(NdEpsBSplinePlatform((3,3)),(0.0..0.5)^2)
P3 = ExtensionFramePlatform(NdCDBSplinePlatform((3,3)),(0.0..0.5)^2)
savefigs = (figname, obj) -> begin
    pgfsave(figname * ".pdf", obj)
    run(`pdf2svg $(figname * ".pdf") $(figname * ".svg")`)
    pgfsave(figname * ".tex", obj);
    pgfsave(figname * ".tikz", obj;include_preamble=false);
    return nothing
end
```

```@meta
DocTestSetup = quote
    using BSplineExtension, DomainSets  
end
```

The package provides a fast solver for multidimensional B-spline extensions such as, e.g.,
```jldoctest Ndframe
julia> P1 = ExtensionFramePlatform(NdBSplinePlatform((3,3)),(0.0..0.5)^2)
ExtensionFramePlatform(ProductPlatform{2}((BSplinePlatform{Float64,3}(), BSplinePlatform{Float64,3}())), 0.0..0.5 x 0.0..0.5)
julia> P2 = ExtensionFramePlatform(NdEpsBSplinePlatform((3,3)),(0.0..0.5)^2)
ExtensionFramePlatform(ProductPlatform{2}((EpsBSplinePlatform{Float64,3}(), EpsBSplinePlatform{Float64,3}())), 0.0..0.5 x 0.0..0.5)
julia> P3 = ExtensionFramePlatform(NdCDBSplinePlatform((3,3)),(0.0..0.5)^2)
ExtensionFramePlatform(ProductPlatform{2}((CDBSplinePlatform{Float64,3}(), CDBSplinePlatform{Float64,3}())), 0.0..0.5 x 0.0..0.5)
```

## Column truncation
In one dimensional B-spline extension approximation the nonzero columns of $$A-AZ'A$$
was independent of $$N$$. It is, however dependent on the size of the boundary of the domain. For two dimensional domains the number of nonzero rows is thus $$\mathcal O(\sqrt{N})$$.

[](
<!-- @example Ndframe
Ns = 6*[1<<k for k in 1:1:4]
ds = 1:4 # hide
PLATFORMs = (NdBSplinePlatform, NdEpsBSplinePlatform, NdCDBSplinePlatform) # hide
crop_tols = 10.0.^(-16.:11.:-5.) # hide
colsizes = Array{Int}(undef, length(PLATFORMs), length(ds), length(Ns), length(crop_tols)) # hide
rowsizes = similar(colsizes) # hide
length(rowsizes)
for (i,PLATFORM) in enumerate(PLATFORMs), (j,d) in enumerate(ds), (k,N) in enumerate(Ns), (l,crop_tol) in enumerate(crop_tols) # hide
    @show d, N
    P = ExtensionFramePlatform(PLATFORM((d,d)), (0.0..0.5)^2); # hide
    M = firstAZstepoperator(P,(N,N);L=(4N,4N)); # hide
    colsizes[i,j,k,l], rowsizes[i,j,k,l]  = truncated_size(BSplineExtensionSolver(M; crop_tol=crop_tol,lazy=true)) # hide
end # hide
A = [];for (i,PLATFORM) in enumerate(PLATFORMs) # hide
    push!(A,@pgf {xmode="log",ymode="log",legend_pos="north west",xlabel="N",title=["BSplinePlatform","EpsBSplinePlatform","CDBSplinePlatform"][i]}) # hide
    for (j,d) in enumerate(ds) # hide
        push!(A, @pgf PlotInc({color=["blue", "red", "brown", "black"][j],mark=["o", "square", "diamond", "x"][j],mark_options="solid"}, Table(Ns.^2, rowsizes[i,j,:,1]))) # hide
        i==1 && push!(A, @pgf LegendEntry("m=$d")) # hide
    end # hide
    push!(A,  @pgf Plot({dashed,color="black"}, Table(Ns.^2, 3Ns))) # hide
    i==1 && push!(A, @pgf LegendEntry(L"\mathcal O(N^{\frac{1}{2}})")) # hide
end # hide
P = @pgf PGFPlotsX.GroupPlot({ymin=0,group_style={group_size="3 by 1",},}, # hide
    A...) # hide
savefigs("2dtruncated_size_1",P)# hide -->
)

[\[.pdf\]](figs/2dtruncated_size_1.pdf), [\[generated .tex\]](figs/2dtruncated_size_1.tex)

![](figs/2dtruncated_size_1.svg)

## Row truncation
Also the number of rows with a norm larger than epsilon is $$\mathcal O(\sqrt{N})$$ since the support of a B-spline is independent of $$N$$. In two dimension, the support is of the dual dictianary (as in `P1` and `P2`) is quite large. This is why the growth in nonzero rows is first $$\mathcal O(N)$$ in the region that is displayed below in the first two panels. Eventually, the  $$\mathcal O(\sqrt{N})$$ will be reached. For the compact dual splines the level off is reached at much lower $$N$$.

[](
 <!-- @example Ndframe
A = [];for (i,PLATFORM) in enumerate(PLATFORMs) # hide
    push!(A,@pgf {xlabel="N",xmode="log",ymode="log",legend_pos="north west",title=["BSplinePlatform","EpsBSplinePlatform","CDBSplinePlatform"][i]}) # hide
    for (j,d) in enumerate(ds) # hide
        opts = @pgf {color=["blue", "red", "brown", "black"][j],mark=["o", "square", "diamond", "x"][j],mark_options="solid"} # hide
        for (l,crop_tol) in enumerate(crop_tols) # hide
            push!(A, @pgf Plot({opts..., (@pgf {solid}, {dashed})[l]...}, Table(Ns.^2, colsizes[i,j,:,l]))) # hide
            i==1 && push!(A, @pgf LegendEntry("m=$d($(crop_tol))")) # hide
        end # hide
    end # hide
    push!(A,  @pgf Plot({dashed,color="black"}, Table(Ns.^2, 500Ns))) # hide
    i==1 && push!(A, @pgf LegendEntry(L"\mathcal O(N^{\frac{1}{2}})")) # hide
end # hide
P = @pgf PGFPlotsX.GroupPlot({ymin=0,group_style={group_size="3 by 1",},}, # hide
    A...) # hide
savefigs("2dtruncated_size_2",P) #hide -->
)


[\[.pdf\]](figs/2dtruncated_size_2.pdf), [\[generated .tex\]](figs/2dtruncated_size_2.tex)

![](figs/2dtruncated_size_2.svg)

# Multidimensional B-spline extension approximation
In this section we use the `BSplineExtensionSolver`](@ref) solver to approximate
the function $$f(x,y)=e^{x*y}$$ on the square $$[0,0.5]^2$$ using a B-spline basis of order
`m` on the interval $$[0,1]^2$$. First we show the convergence results, then the time complexity
of the approximation algorithm.
## Convergence
In the figure below, which shows the uniform error of approximating a 2 dimensional analytic function (details are in the introduction of this section),
 we see that convergence is algebraic.

[](
<!-- @example Ndframe
Ns1 = 6*[1<<k for k in 1:1:4]
Ns2 = 6*[1<<k for k in 1:1:4]
ds = 1:4 # hide
errors = Array{Float64}(undef, length(PLATFORMs), length(ds), length(Ns1)) # hide
timings = Array{Float64}(undef, length(PLATFORMs), length(ds), length(Ns2)) # hide
for (i,PLATFORM) in enumerate(PLATFORMs), (j,d) in enumerate(ds), (k,N) in enumerate(Ns1) #hide
    P = ExtensionFramePlatform(PLATFORM((d,d)), (0.0..0.5)^2) #hide
    F,_ = @timed Fun(f, P, (N,N);L=(4N,4N), REG=BSplineExtension.BSplineExtensionSolver, crop=true, crop_tol=1e-10) #hide
    errors[i,j,k] = abserror(f, F) #hide
end #hide
for (i,PLATFORM) in enumerate(PLATFORMs), (j,d) in enumerate(ds), (k,N) in enumerate(Ns2) #hide
    P = ExtensionFramePlatform(PLATFORM((d,d)), (0.0..0.5)^2) #hide
    _,timings[i,j,k],_ = @timed Fun(f, P, (N,N);L=(4N,4N), REG=BSplineExtension.BSplineExtensionSolver, crop=true, crop_tol=1e-10) #hide
end # hide
A = [];for (i,PLATFORM) in enumerate(PLATFORMs) # hide
    push!(A,@pgf {xmode="log",ymode="log",xlabel="N",legend_pos="north west",title=["BSplinePlatform","EpsBSplinePlatform","CDBSplinePlatform"][i]}) # hide
    for (j,d) in enumerate(ds) # hide
        opts = @pgf {color=["blue", "red", "brown", "black"][j],mark=["o", "square", "diamond", "x"][j],mark_options="solid"} # hide
        push!(A, @pgf Plot(opts, Table(Ns1.^2, errors[i,j,:]))) # hide
        i==1 && push!(A, @pgf LegendEntry("m=$d")) # hide
    end # hide
end # hide
P = @pgf PGFPlotsX.GroupPlot({ymin=0,group_style={group_size="3 by 1",},}, # hide
    A...) # hide
savefigs("2derrors",P) # hide  -->
)


[\[.pdf\]](figs/2derrors.pdf), [\[generated .tex\]](figs/2derrors.tex)

![](figs/2derrors.svg)

## Timings
For all platforms and for high $$N$$, the most costly part is the first step of the AZ algorithm, i.e., solving $$(A-AZ'A)x=(I-AZ')b$$.
The system is of size $$\mathcal O(n\times n) = \mathcal O(\sqrt{N}\times\sqrt{N})$$, where $$N$$ is the total number of
degrees of freedom and $$n$$ is the number of dof in one dimension. Therefore, solving this system and by consequence the approting a function costs
$$\mathcal O(n^3)=\mathcal O(N^{\frac{3}{2}})$$.

For low $$N$$ we saw that the first two platforms are not in the regime where the size of the system of $$\mathcal O(\sqrt{N}\times\sqrt{N})$$, they are $$\mathcal O({N}\times\sqrt{N})$$. We expect to see a $$\mathcal O(N^{2})$$ complexity. This complexity
is for small $$N$$ overshadowed by the $$\mathcal O({N^{\tfrac{3}{2}}})$$ cost of creating the system of
the first AZ step.

[](
<!-- example Ndframe
A = [];for (i,PLATFORM) in enumerate(PLATFORMs) # hide
    push!(A,@pgf {xmode="log",ymode="log",xlabel="N",legend_pos="north west",title=["BSplinePlatform","EpsBSplinePlatform","CDBSplinePlatform"][i]}) # hide
    for (j,d) in enumerate(ds) # hide
        opts = @pgf {color=["blue", "red", "brown", "black"][j],mark=["o", "square", "diamond", "x"][j],mark_options="solid"} # hide
        push!(A, @pgf Plot(opts, Table(Ns2.^2, timings[i,j,:]))) # hide
        i==1 && push!(A, @pgf LegendEntry("m=$d")) # hide
    end # hide
    push!(A, @pgf Plot({color="black",dashed},Table(Ns2.^2,5e-5Ns2.^2(1.5)))) # hide
    i==1 && push!(A, LegendEntry(L"\mathcal O(N^\frac{3}{2})")) # hide
end # hide
P = @pgf PGFPlotsX.GroupPlot({ymin=0,group_style={group_size="3 by 1",},}, # hide
    A...) # hide
savefigs("2dtimings",P) # hide -->
)

[\[.pdf\]](figs/2dtimings.pdf), [\[generated .tex\]](figs/2dtimings.tex)

![](figs/2dtimings.svg)
