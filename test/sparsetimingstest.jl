using BSplineExtension, StaticArrays, LinearAlgebra, SparseArrays
using BSplineExtension.AZSparse: sparsemixedgramcomplement, sparseRAE
using BSplineExtension.BSplineExtensionSolvers: nonzero_rows
D = .4*ball() + SVector(.5,.5,.5)
Pbasis = NdCDBSplinePlatform((2,2,2))
P = ExtensionFramePlatform(Pbasis, D)
m = (2,2,2)
ns = 10:10:200

n = 5
N = (n,n,n)
dict1 = dictionary(P,N)
g = oversampling_grid(P,N;L=m.*N)
μ = discretemeasure(g)
dict2 = dualdictionary(P,N,μ)

ImZA, t, _= @timed sparsemixedgramcomplement(dict2, dict1, μ)
col_indices, _ = @timed findall(reshape(nonzero_rows(ImZA),size(dict1)))
A1, t, _ = @timed sparseRAE(dict1, g, col_indices)
A2, _ = @timed ImZA[LinearIndices(size(dict1))[col_indices],:]
A, t, _ = @timed A1*A2
QRfact,t,_ = @timed qr(A);
Base.summarysize(QRfact.Q) + Base.summarysize(QRfact.R) + Base.summarysize(QRfact.pcol) + Base.summarysize(QRfact.prow)
Base.summarysize(QRfact)


timings_stepA2 = zeros(length(ns))
timings_stepA1 = similar(timings_stepA2)
timings_stepA1A2 = similar(timings_stepA2)
timings_stepQR = similar(timings_stepA2)
nzeroels = Vector{Int}(undef,length(ns))
nnzels = Vector{Int}(undef,length(ns))
sQ = Vector{Int}(undef,length(ns))
sR = Vector{Int}(undef,length(ns))
sT = Vector{Int}(undef,length(ns))
sProw = Vector{Int}(undef,length(ns))
sPcol = Vector{Int}(undef,length(ns))
for (i,n) in enumerate(ns)
    N = (n,n,n)
    @show N
    dict1 = dictionary(P,N)
    g = oversampling_grid(P,N;L=m.*N)
    μ = discretemeasure(g)
    dict2 = dualdictionary(P,N,μ)

    ImZA, t, _= @timed sparsemixedgramcomplement(dict2, dict1, μ)
    @show timings_stepA2[i] = t
    col_indices = findall(reshape(nonzero_rows(ImZA),size(dict1)))
    A1, t, _ = @timed sparseRAE(dict1, g, col_indices)
    @show timings_stepA1[i] = t
    A2 = ImZA[LinearIndices(size(dict1))[col_indices],:]
    A, t, _ = @timed A1*A2
    @show timings_stepA1A2[i] = t
    @show nzeroels[i] = count(abs.(A.nzval) .< 1e-12)
    @show nnzels[i] = nnz(A)
    QRfact,t,_ = @timed qr(A)
    @show timings_stepQR[i] = t
    sQ[i], sT[i], sR[i], sPcol[i],sProw[i] = map(x->Base.summarysize(getfield(QRfact,x)), fieldnames(typeof(QRfact)))
end
@show ns
@show timings_stepA2
@show timings_stepA1
@show timings_stepA1A2
@show timings_stepQR
@show nzeroels
@show nnzels
@show sQ
@show sR
@show sT
@show sPcol
@show sProw

# ns = 10:10:100
# timings_stepA2 = [0.230505, 0.91782, 2.08745, 3.93506, 6.33685, 9.32926, 12.7722, 16.185, 20.8174, 25.6926]
# timings_stepA1 = [0.0032484, 0.100566, 0.114128, 0.135143, 0.095328, 0.135276, 0.202571, 0.247372, 0.339048, 0.408019]
# timings_stepA1A2 = [0.0116177, 0.0580797, 0.182542, 0.185742, 0.330227, 0.556001, 0.69891, 0.89139, 1.20048, 1.3486]
# timings_stepQR = [0.0873553, 2.07531, 5.45822, 12.7753, 34.9328, 62.6667, 108.627, 186.835, 396.584, 441.746]
# nzeroels = [0, 288, 3408, 4896, 9648, 13920, 17160, 24624, 31176, 41640]
# nnzels = [242822, 1197634, 2798866, 5418074, 8798074, 12436402, 17424026, 22703930, 28742482, 35870426]
# sQ = [15360064, 265807272, 707821936, 1504745896, 3238940768, 4945652176, 7157216392, 10436664312, 16911629752, 18862633584]
# sR = [2558080, 33734488, 92452552, 198993320, 421433240, 637232040, 926261688, 1332687720, 2156718504, 2388866632]
# sT = [4512, 86440, 309840, 672696, 1100560, 1646512, 2395736, 3149208, 4040024, 5179168]
# sPcol = [4704, 19128, 42552, 76728, 120696, 172344, 235704, 308280, 389304, 482136]
# sProw = [16912, 136656, 462256, 1096560, 2142128, 3702288, 5882288, 8783376, 12499760, 17149072]
#
#
# using PGFPlotsX, LaTeXStrings
# @pgf Axis({xmode="log",ymode="log"},
#     Plot({only_marks},Table(nnzels,ns.^3))
#
#     )
# @pgf Axis({legend_pos="north west",legend_cell_align="left",xmode="log",ymode="log"},
#     Plot({only_marks,color="blue"},Table(nnzels,timings_stepA2)),
#     LegendEntry(L"I-AZ^*"),
#     Plot({only_marks,color="red"},Table(nnzels,timings_stepA1)),
#     LegendEntry(L"A"),
#     Plot({only_marks,color="green"},Table(nnzels,timings_stepA1A2)),
#     LegendEntry(L"\times"),
#     Plot({only_marks,color="black"},Table(nnzels,timings_stepQR)),
#     LegendEntry(L"QR"),
#     )
# P = @pgf Axis({xlabel="N",legend_pos="north west",legend_cell_align="left",xmode="log",ymode="log"},
#     Plot({only_marks,color="blue"},Table(nnzels,sQ)),
#     LegendEntry("factors (Q)"),
#     Plot({only_marks,color="red"},Table(nnzels,sT)),
#     LegendEntry(L"\tau"),
#     Plot({only_marks,color="green"},Table(nnzels,sR)),
#     LegendEntry(L"R"),
#     Plot({only_marks,color="black"},Table(nnzels,sPcol)),
#     LegendEntry("cpiv"),
#     Plot({only_marks,color="brown"},Table(nnzels,sProw)),
#     LegendEntry("rpivinv"),
#     Plot({style="black,dashed"},Table(nnzels,nnzels)),
#     LegendEntry(L"\mathcal O(N)"),
# )
# using DocumentPGFPlots
# DocumentPGFPlots.savefigs(joinpath(pwd(),"memorycheck",),P)
#
