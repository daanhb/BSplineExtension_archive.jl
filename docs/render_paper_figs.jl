using BSplineExtension, DomainSets, PGFPlotsX, DocumentPGFPlots, LaTeXStrings
@info "Start rendering paper figures"
imgdir = normpath(joinpath(@__DIR__(),"src","man","figs"))
all=false

Ns = 6*[1<<k for k in 1:1:4]
ds = 1:4 # hide
PLATFORMs = (NdBSplinePlatform, NdCDBSplinePlatform) # hide
colsizes = Array{Int}(undef, length(PLATFORMs), length(ds), length(Ns)) # hide
rowsizes = similar(colsizes) # hide
length(rowsizes)
for (i,PLATFORM) in enumerate(PLATFORMs), (j,d) in enumerate(ds), (k,N) in enumerate(Ns) # hide
    @show d, N
    P = ExtensionFramePlatform(PLATFORM((d,d)), (0.0..0.5)^2); # hide
    M = firstAZstepoperator(P,(N,N)); # hide
    colsizes[i,j,k], rowsizes[i,j,k]  = truncated_size(BSplineExtensionSolver(M; crop_tol=1e-5,lazy=true)) # hide
end # hide
A = [];
    push!(A,@pgf {xmode="log",ymode="log",legend_pos="north west",xlabel="N",ymin=1e1,ymax=3e4}) # hide
    for (j,d) in enumerate(ds) # hide
        push!(A, @pgf PlotInc({color=["blue", "red", "brown", "black"][j],mark=["o", "square", "diamond", "x"][j],mark_options="solid"}, Table(Ns.^2, rowsizes[1,j,:]))) # hide
        push!(A, @pgf LegendEntry("m=$d")) # hide
    end # hide
    push!(A,  @pgf Plot({dashed,color="black"}, Table(Ns.^2, 3Ns))) # hide
    push!(A, @pgf LegendEntry(L"\mathcal O(N^{\frac{1}{2}})")) # hide
    push!(A,@pgf {xmode="log",ymode="log",xlabel="N",ymin=1e1,ymax=3e4}) # hide
    for (j,d) in enumerate(ds) # hide
        push!(A, @pgf PlotInc({color=["blue", "red", "brown", "black"][j],mark=["o", "square", "diamond", "x"][j],mark_options="solid"}, Table(Ns.^2, colsizes[1,j,:]))) # hide
    end # hide
    push!(A, @pgf Plot({dashed,color="black"}, Table(Ns.^2, 3Ns))) # hide
    push!(A, @pgf {xmode="log",ymode="log",xlabel="N",ymin=1e1,ymax=3e4}) # hide
    for (j,d) in enumerate(ds) # hide
        push!(A, @pgf PlotInc({color=["blue", "red", "brown", "black"][j],mark=["o", "square", "diamond", "x"][j],mark_options="solid"}, Table(Ns.^2, colsizes[2,j,:]))) # hide
    end # hide
    push!(A,  @pgf Plot({dashed,color="black"}, Table(Ns.^2, 3Ns))) # hide

P = @pgf PGFPlotsX.GroupPlot({ymin=0,group_style={group_size="3 by 1",},}, # hide
    A...) # hide
DocumentPGFPlots.savefigs(joinpath(imgdir,"truncatedsizepaper"),P;all=all)





Ns1 = 6*[1<<k for k in 1:1:4]
Ns2 = 6*[1<<k for k in 1:1:4]
verbose=true
ds = 1:4 # hide
errors = Array{Float64}(undef, 3, length(ds), length(Ns1)) # hide
timings = Array{Float64}(undef, 3, length(ds), length(Ns2)) # hide
f = (x,y)->exp(x*y)
i = 1
    for (j,d) in enumerate(ds), (k,N) in enumerate(Ns1) #hide
        P = ExtensionFramePlatform(NdEpsBSplinePlatform((d,d)), (0.0..0.5)^2) #hide
        F,_ = @timed Fun(f, P, (N,N);REG=BSplineExtension.BSplineExtensionSolver, crop=true, crop_tol=1e-10, verbose=verbose) #hide
        _,timings[i,j,k],_ = @timed Fun(f, P, (N,N);REG=BSplineExtension.BSplineExtensionSolver, crop=true, crop_tol=1e-10) #hide
        errors[i,j,k] = abserror(f, F) #hide
    end #hide
    i = 2
    for (j,d) in enumerate(ds), (k,N) in enumerate(Ns1) #hide
        P = ExtensionFramePlatform(NdCDBSplinePlatform((d,d)), (0.0..0.5)^2) #hide
        F,_ = @timed Fun(f, P, (N,N);REG=BSplineExtension.BSplineExtensionSolver, crop=true, crop_tol=1e-10, verbose=verbose) #hide
        _,timings[i,j,k],_ = @timed Fun(f, P, (N,N);REG=BSplineExtension.BSplineExtensionSolver, crop=true, crop_tol=1e-10) #hide
        errors[i,j,k] = abserror(f, F) #hide
    end #hide
    i = 3
    for (j,d) in enumerate(ds), (k,N) in enumerate(Ns1) #hide
        P = ExtensionFramePlatform(NdCDBSplinePlatform((d,d)), (0.0..0.5)^2) #hide
        F,_ = @timed Fun(f, P, (N,N);REG=BSplineExtension.BSplineExtensionSolver, crop=true, crop_tol=1e-10, sparse=true,verbose=verbose) #hide
        _,timings[i,j,k],_ =  @timed Fun(f, P, (N,N);REG=BSplineExtension.BSplineExtensionSolver, crop=true, crop_tol=1e-10, sparse=true) #hide
        errors[i,j,k] = abserror(f, F) #hide
    end #hide

A = [];for i in 1:3 # hide
    push!(A,@pgf {xmode="log",ymode="log",xlabel="N",legend_pos="north west"}) # hide
    for (j,d) in enumerate(ds) # hide
        opts = @pgf {color=["blue", "red", "brown", "black"][j],mark=["o", "square", "diamond", "x"][j],mark_options="solid"} # hide
        push!(A, @pgf Plot(opts, Table(Ns1.^2, errors[i,j,:]))) # hide
        i==1 && push!(A, @pgf LegendEntry("m=$d")) # hide
    end # hide
end # hide
P = @pgf PGFPlotsX.GroupPlot({ymin=0,group_style={group_size="3 by 1",},}, # hide
    A...) # hide
DocumentPGFPlots.savefigs(joinpath(imgdir,"2derrorpaper"),P;all=all) #

A = [];for i in 1:3  # hide
    push!(A,@pgf {xmode="log",ymode="log",xlabel="N",legend_pos="north west",ymin=1e-2,ymax=1e2}) # hide
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
DocumentPGFPlots.savefigs(joinpath(imgdir,"2dtimingspaper"),P;all=all)




f = (x,y,z)->exp(x*y*y)
Ns1 = round.(Int,exp.(LinRange(log.((1000 ,100000))...,4)).^(.33))
Ns2 = round.(Int,exp.(LinRange(log.((1000 ,100000))...,4)).^(.33))
ds = 1:4 # hide
errors = Array{Float64}(undef, 2, length(ds), length(Ns1)) # hide
timings = Array{Float64}(undef, 2, length(ds), length(Ns2)) # hide
i = 1
    for (j,d) in enumerate(ds), (k,N) in enumerate(Ns1) #hide
        @show d, N
        P = ExtensionFramePlatform(NdCDBSplinePlatform((d,d,d)), (0.0..0.5)^3) #hide
        F,_ = @timed Fun(f, P, (N,N,N);REG=BSplineExtension.BSplineExtensionSolver, crop=true, crop_tol=1e-10, verbose=verbose) #hide
        _,timings[i,j,k],_ = @timed Fun(f, P, (N,N,N);REG=BSplineExtension.BSplineExtensionSolver, crop=true, crop_tol=1e-10) #hide
        errors[i,j,k] = abserror(f, F) #hide
    end #hide
i = 2
    for (j,d) in enumerate(ds), (k,N) in enumerate(Ns1) #hide
        @show d, N
        P = ExtensionFramePlatform(NdCDBSplinePlatform((d,d,d)), (0.0..0.5)^3) #hide
        F,_ = @timed Fun(f, P, (N,N,N);REG=BSplineExtension.BSplineExtensionSolver, crop=true, crop_tol=1e-10, sparse=true, verbose=verbose) #hide
        errors[i,j,k] = abserror(f, F) #hide
        _,timings[i,j,k],_ = @timed Fun(f, P, (N,N,N);REG=BSplineExtension.BSplineExtensionSolver, crop=true, sparse=true, crop_tol=1e-10) #hide
    end #hide

A = [];for i in 1:2 # hide
    push!(A,@pgf {xmode="log",ymode="log",xlabel="N",legend_pos="north west"}) # hide
    for (j,d) in enumerate(ds) # hide
        opts = @pgf {color=["blue", "red", "brown", "black"][j],mark=["o", "square", "diamond", "x"][j],mark_options="solid"} # hide
        push!(A, @pgf Plot(opts, Table(Ns1.^2, errors[i,j,:]))) # hide
        i==1 && push!(A, @pgf LegendEntry("m=$d")) # hide
    end # hide
end # hide
P = @pgf PGFPlotsX.GroupPlot({ymin=0,group_style={group_size="2 by 1",},}, # hide
    A...) # hide
DocumentPGFPlots.savefigs(joinpath(imgdir,"3derrorpaper"),P;all=all) #

A = [];for i in 1:2  # hide
    push!(A,@pgf {xmode="log",ymode="log",xlabel="N",legend_pos="north west"}) # hide
    for (j,d) in enumerate(ds) # hide
        opts = @pgf {color=["blue", "red", "brown", "black"][j],mark=["o", "square", "diamond", "x"][j],mark_options="solid"} # hide
        push!(A, @pgf Plot(opts, Table(Ns2.^3, timings[i,j,:]))) # hide
        i==1 && push!(A, @pgf LegendEntry("m=$d")) # hide
    end # hide
    push!(A, @pgf Plot({color="black",dashed},Table(Ns2.^3,5e-5Ns2.^6))) # hide
    i==1 && push!(A, LegendEntry(L"\mathcal O(N^2)")) # hide
end # hide
P = @pgf PGFPlotsX.GroupPlot({ymin=0,group_style={group_size="2 by 1",},}, # hide
    A...) # hide
DocumentPGFPlots.savefigs(joinpath(imgdir,"3dtimingspaper"),P;all=all)


DocumentPGFPlots.savefigs(joinpath(imgdir,"truncatedsizepaper"),P;render=false)
DocumentPGFPlots.savefigs(joinpath(imgdir,"2derrorpaper"),P;render=false)
DocumentPGFPlots.savefigs(joinpath(imgdir,"2dtimingspaper"),P;render=false)
DocumentPGFPlots.savefigs(joinpath(imgdir,"3derrorpaper"),P;render=false) #
DocumentPGFPlots.savefigs(joinpath(imgdir,"3dtimingspaper"),P;render=false)
