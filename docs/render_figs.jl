using BSplineExtension, DomainSets, PGFPlotsX, DocumentPGFPlots, LaTeXStrings
@info "Start rendering figures"
imgdir = normpath(joinpath(@__DIR__(),"src","man","figs"))

P1 = BSplinePlatform()
P2 = EpsBSplinePlatform()
P3 = CDBSplinePlatform()
f = x->exp(cos(10pi*x));N=60;
F1 = Fun(f, P1, N)
F2 = Fun(f, P2, N)
F3 = Fun(f, P3, N)
P = @pgf GroupPlot({group_style = {group_size="2 by 1",},},
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
        PlotInc(f, F3;n=600))
DocumentPGFPlots.savefigs(joinpath(imgdir,"basis"), P) # hide

Ns = 10*[1<<k for k in 0:11] # hide
err = zeros(3, length(Ns)) # hide
for (i,P) in enumerate((P1,P2,P3))
    for (j,N) in enumerate(Ns)
        F = Fun(f, P, N)
        err[i,j] = abserror(f, F)
    end
end
P = @pgf Axis({ymode="log",xmode="log"},
        PlotInc(Table(Ns,err[1,:])),
        LegendEntry("P1"),
        PlotInc(Table(Ns,err[2,:])),
        LegendEntry("P2"),
        PlotInc(Table(Ns,err[3,:])),
        LegendEntry("P3"),
        )
DocumentPGFPlots.savefigs(joinpath(imgdir,"convergence_basis"), P) # hide


Ns = 20:20:300 # hide
ds = 1:4 # hide
PLATFORMs = (BSplinePlatform, EpsBSplinePlatform, CDBSplinePlatform) # hide
crop_tols = 10.0.^(-16.:6.:-10.) # hide
colsizes = Array{Int}(undef, length(PLATFORMs), length(ds), length(Ns), length(crop_tols)) # hide
rowsizes = similar(colsizes) # hide
for (i,PLATFORM) in enumerate(PLATFORMs), (j,d) in enumerate(ds), (k,N) in enumerate(Ns), (l,crop_tol) in enumerate(crop_tols) # hide
    P = ExtensionFramePlatform(PLATFORM(d), 0.0..0.5); # hide
    plunge = plungeoperator(P,N;); A = AZ_A(P,N;); Zt = AZ_Zt(P,N); # hide
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
P = @pgf PGFPlotsX.GroupPlot({ymin=0,group_style={group_size="3 by 1",},}, # hide
    A...) # hide
DocumentPGFPlots.savefigs(joinpath(imgdir,"truncated_size_1"), P) # hide

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
P = @pgf PGFPlotsX.GroupPlot({ymin=0,group_style={group_size="3 by 1",},}, # hide
    A...) # hide
DocumentPGFPlots.savefigs(joinpath(imgdir,"truncated_size_2"), P) # hide

PLATFORMs = (EpsBSplinePlatform, BSplinePlatform, CDBSplinePlatform) # hide
Ns1 = [1<<k for k in 4:10] # hide
Ns2 = [1<<k for k in 9:16] # hide
ds = 1:4 # hide
errors = Array{Float64}(undef, length(PLATFORMs), length(ds), length(Ns1)) # hide
timings = Array{Float64}(undef, length(PLATFORMs), length(ds), length(Ns2)) # hide
for (i,PLATFORM) in enumerate(PLATFORMs), (j,d) in enumerate(ds), (k,N) in enumerate(Ns1) #hide
    P = ExtensionFramePlatform(PLATFORM(d), 0.0..0.5) #hide
    F,_ = @timed Fun(exp, P, N;REG=BSplineExtension.BSplineExtensionSolver, crop=true, crop_tol=1e-10) #hide
    errors[i,j,k] = abserror(exp, F) #hide
end #hide
for (i,PLATFORM) in enumerate(PLATFORMs), (j,d) in enumerate(ds), (k,N) in enumerate(Ns2) #hide
    P = ExtensionFramePlatform(PLATFORM(d), 0.0..0.5); #hide
    _,timings[i,j,k],_ = @timed Fun(exp, P, N;REG=BSplineExtension.BSplineExtensionSolver,  crop=true, crop_tol=1e-10) #hide
end # hide
A = [];for (i,PLATFORM) in enumerate(PLATFORMs) # hide
    push!(A,@pgf {xmode="log",ymode="log",xlabel="N",legend_pos="north west",title=["BSplinePlatform","EpsBSplinePlatform","CDBSplinePlatform"][i]}) # hide
    for (j,d) in enumerate(ds) # hide
        opts = @pgf {color=["blue", "red", "brown", "black"][j],mark=["o", "square", "diamond", "x"][j],mark_options="solid"} # hide
        push!(A, @pgf Plot(opts, Table(Ns1, errors[i,j,:]))) # hide
        i==1 && push!(A, @pgf LegendEntry("m=$d")) # hide
    end # hide
end # hide
P = @pgf PGFPlotsX.GroupPlot({ymin=0,group_style={group_size="3 by 1",},}, # hide
    A...) # hide
DocumentPGFPlots.savefigs(joinpath(imgdir,"1derrors"), P) # hide


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
P = @pgf PGFPlotsX.GroupPlot({ymin=0,group_style={group_size="3 by 1",},}, # hide
    A...) # hide
DocumentPGFPlots.savefigs(joinpath(imgdir,"1dtimings"), P) # hide


f = (x,y)->exp(x*y)
# Ns = 6*[1<<k for k in 1:1:4]
# ds = 1:4 # hide
# PLATFORMs = (NdBSplinePlatform, NdEpsBSplinePlatform, NdCDBSplinePlatform) # hide
# crop_tols = 10.0.^(-16.:11.:-5.) # hide
# colsizes = Array{Int}(undef, length(PLATFORMs), length(ds), length(Ns), length(crop_tols)) # hide
# rowsizes = similar(colsizes) # hide
# length(rowsizes)
# for (i,PLATFORM) in enumerate(PLATFORMs), (j,d) in enumerate(ds), (k,N) in enumerate(Ns), (l,crop_tol) in enumerate(crop_tols) # hide
#     @show d, N
#     P = ExtensionFramePlatform(PLATFORM((d,d)), (0.0..0.5)^2); # hide
#     M = firstAZstepoperator(P,(N,N)); # hide
#     colsizes[i,j,k,l], rowsizes[i,j,k,l]  = truncated_size(BSplineExtensionSolver(M; crop_tol=crop_tol,lazy=true)) # hide
# end # hide
# A = [];for (i,PLATFORM) in enumerate(PLATFORMs) # hide
#     push!(A,@pgf {xmode="log",ymode="log",legend_pos="north west",xlabel="N",title=["BSplinePlatform","EpsBSplinePlatform","CDBSplinePlatform"][i]}) # hide
#     for (j,d) in enumerate(ds) # hide
#         push!(A, @pgf PlotInc({color=["blue", "red", "brown", "black"][j],mark=["o", "square", "diamond", "x"][j],mark_options="solid"}, Table(Ns.^2, rowsizes[i,j,:,1]))) # hide
#         i==1 && push!(A, @pgf LegendEntry("m=$d")) # hide
#     end # hide
#     push!(A,  @pgf Plot({dashed,color="black"}, Table(Ns.^2, 3Ns))) # hide
#     i==1 && push!(A, @pgf LegendEntry(L"\mathcal O(N^{\frac{1}{2}})")) # hide
# end # hide
# P = @pgf PGFPlotsX.GroupPlot({ymin=0,group_style={group_size="3 by 1",},}, # hide
#     A...) # hide
# DocumentPGFPlots.savefigs(joinpath(imddir,"2dtruncated_size_1"),P)



# A = [];for (i,PLATFORM) in enumerate(PLATFORMs) # hide
#     push!(A,@pgf {xlabel="N",xmode="log",ymode="log",legend_pos="north west",title=["BSplinePlatform","EpsBSplinePlatform","CDBSplinePlatform"][i]}) # hide
#     for (j,d) in enumerate(ds) # hide
#         opts = @pgf {color=["blue", "red", "brown", "black"][j],mark=["o", "square", "diamond", "x"][j],mark_options="solid"} # hide
#         for (l,crop_tol) in enumerate(crop_tols) # hide
#             push!(A, @pgf Plot({opts..., (@pgf {solid}, {dashed})[l]...}, Table(Ns.^2, colsizes[i,j,:,l]))) # hide
#             i==1 && push!(A, @pgf LegendEntry("m=$d($(crop_tol))")) # hide
#         end # hide
#     end # hide
#     push!(A,  @pgf Plot({dashed,color="black"}, Table(Ns.^2, 500Ns))) # hide
#     i==1 && push!(A, @pgf LegendEntry(L"\mathcal O(N^{\frac{1}{2}})")) # hide
# end # hide
# P = @pgf PGFPlotsX.GroupPlot({ymin=0,group_style={group_size="3 by 1",},}, # hide
#     A...) # hide
# DocumentPGFPlots.savefigs(joinpath(imgdir,"2dtruncated_size_2"),P)



# Ns1 = 6*[1<<k for k in 1:1:4]
# Ns2 = 6*[1<<k for k in 1:1:4]
# ds = 1:4 # hide
# errors = Array{Float64}(undef, length(PLATFORMs), length(ds), length(Ns1)) # hide
# timings = Array{Float64}(undef, length(PLATFORMs), length(ds), length(Ns2)) # hide
# for (i,PLATFORM) in enumerate(PLATFORMs), (j,d) in enumerate(ds), (k,N) in enumerate(Ns1) #hide
#     P = ExtensionFramePlatform(PLATFORM((d,d)), (0.0..0.5)^2) #hide
#     F,_ = @timed Fun(f, P, (N,N);REG=BSplineExtension.BSplineExtensionSolver, crop=true, crop_tol=1e-10) #hide
#     errors[i,j,k] = abserror(f, F) #hide
# end #hide
# for (i,PLATFORM) in enumerate(PLATFORMs), (j,d) in enumerate(ds), (k,N) in enumerate(Ns2) #hide
#     P = ExtensionFramePlatform(PLATFORM((d,d)), (0.0..0.5)^2) #hide
#     _,timings[i,j,k],_ = @timed Fun(f, P, (N,N);REG=BSplineExtension.BSplineExtensionSolver, crop=true, crop_tol=1e-10) #hide
# end # hide
# A = [];for (i,PLATFORM) in enumerate(PLATFORMs) # hide
#     push!(A,@pgf {xmode="log",ymode="log",xlabel="N",legend_pos="north west",title=["BSplinePlatform","EpsBSplinePlatform","CDBSplinePlatform"][i]}) # hide
#     for (j,d) in enumerate(ds) # hide
#         opts = @pgf {color=["blue", "red", "brown", "black"][j],mark=["o", "square", "diamond", "x"][j],mark_options="solid"} # hide
#         push!(A, @pgf Plot(opts, Table(Ns1.^2, errors[i,j,:]))) # hide
#         i==1 && push!(A, @pgf LegendEntry("m=$d")) # hide
#     end # hide
# end # hide
# P = @pgf PGFPlotsX.GroupPlot({ymin=0,group_style={group_size="3 by 1",},}, # hide
#     A...) # hide
# DocumentPGFPlots.savefigs(joinpath(imddir,"2derrors"),P) #



# A = [];for (i,PLATFORM) in enumerate(PLATFORMs) # hide
#     push!(A,@pgf {xmode="log",ymode="log",xlabel="N",legend_pos="north west",title=["BSplinePlatform","EpsBSplinePlatform","CDBSplinePlatform"][i]}) # hide
#     for (j,d) in enumerate(ds) # hide
#         opts = @pgf {color=["blue", "red", "brown", "black"][j],mark=["o", "square", "diamond", "x"][j],mark_options="solid"} # hide
#         push!(A, @pgf Plot(opts, Table(Ns2.^2, timings[i,j,:]))) # hide
#         i==1 && push!(A, @pgf LegendEntry("m=$d")) # hide
#     end # hide
#     push!(A, @pgf Plot({color="black",dashed},Table(Ns2.^2,5e-5Ns2.^2(1.5)))) # hide
#     i==1 && push!(A, LegendEntry(L"\mathcal O(N^\frac{3}{2})")) # hide
# end # hide
# P = @pgf PGFPlotsX.GroupPlot({ymin=0,group_style={group_size="3 by 1",},}, # hide
#     A...) # hide
# DocumentPGFPlots.savefigs(joinpath(imgdir,"2dtimings"),P)
@info "End rendering figures"
