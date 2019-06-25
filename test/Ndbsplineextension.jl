using Pkg;Pkg.activate(localprojectdir())
using BSplineExtension, Test, LinearAlgebra

using PGFPlotsX, LaTeXStrings
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
    i==1 && push!(A, @pgf LegendEntry(L"\mathcal O(N^{\frac{1}{2}})"))
end # hide
P = @pgf PGFPlotsX.GroupPlot({ymin=0,group_style={group_size="3 by 1",},}, # hide
    A...) # hide
savefigs(joinpath(homedir(), "julia","BSplineExtension.jl","docs","src","man","figs","2dtruncated_size_1"),P)
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
    i==1 && push!(A, @pgf LegendEntry(L"\mathcal O(N^{\frac{1}{2}})"))
end # hide
P = @pgf PGFPlotsX.GroupPlot({ymin=0,group_style={group_size="3 by 1",},}, # hide
    A...) # hide
savefigs(joinpath(homedir(), "julia","BSplineExtension.jl","docs","src","man","figs","2dtruncated_size_2"),P)


using PGFPlotsX, LaTeXStrings
savefigs = (figname, obj) -> begin
    pgfsave(figname * ".pdf", obj)
    run(`pdf2svg $(figname * ".pdf") $(figname * ".svg")`)
    pgfsave(figname * ".tex", obj);
    return nothing
end
f = (x,y)->exp(x*y)
PLATFORMs = (NdEpsBSplinePlatform, NdBSplinePlatform, NdCDBSplinePlatform) # hide
# Ns1 = [1<<k for k in 4:10] # hide
# Ns2 = [1<<k for k in 9:16] # hide
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
savefigs(joinpath(homedir(), "julia","BSplineExtension.jl","docs","src","man","figs","2derrors"),P)
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
savefigs(joinpath(homedir(), "julia","BSplineExtension.jl","docs","src","man","figs","2dtimings"),P)
