

# using BSplineExtension, StaticArrays, LinearAlgebra, SparseArrays
# using BSplineExtension.AZSparse: sparsemixedgramcomplement, sparseRAE
# using BSplineExtension.BSplineExtensionSolvers: nonzero_rows
# D = .4*disk() + SVector(.5,.5)
# Pbasis = NdCDBSplinePlatform((2,2))
# P = ExtensionFramePlatform(Pbasis, D)
# N = 20,20
# m = (2,2)
# ns = 100:100:4000
#
# n = 100
# N = (n,n)
# dict1 = dictionary(P,N)
# g = oversampling_grid(P,N;L=m.*N)
# μ = discretemeasure(g)
# dict2 = dualdictionary(P,N,μ)
#
# ImZA, t, _= @timed sparsemixedgramcomplement(dict2, dict1, μ)
# col_indices, _ = @timed findall(reshape(nonzero_rows(ImZA),size(dict1)))
# A1, t, _ = @timed sparseRAE(dict1, g, col_indices)
# A2, _ = @timed ImZA[LinearIndices(size(dict1))[col_indices],:]
# A, t, _ = @timed A1*A2
# QRfact,t,_ = @timed qr(A);
# Base.summarysize(QRfact.Q) + Base.summarysize(QRfact.R) + Base.summarysize(QRfact.pcol) + Base.summarysize(QRfact.prow)
# Base.summarysize(QRfact)
#
#
# timings_stepA2 = zeros(length(ns))
# timings_stepA1 = similar(timings_stepA2)
# timings_stepA1A2 = similar(timings_stepA2)
# timings_stepQR = similar(timings_stepA2)
# nzeroels = Vector{Int}(undef,length(ns))
# nnzels = Vector{Int}(undef,length(ns))
# sQ = Vector{Int}(undef,length(ns))
# sR = Vector{Int}(undef,length(ns))
# sT = Vector{Int}(undef,length(ns))
# sProw = Vector{Int}(undef,length(ns))
# sPcol = Vector{Int}(undef,length(ns))
# for (i,n) in enumerate(ns)
#     N = (n,n)
#     @show N
#     dict1 = dictionary(P,N)
#     g = oversampling_grid(P,N;L=m.*N)
#     μ = discretemeasure(g)
#     dict2 = dualdictionary(P,N,μ)
#
#     ImZA, t, _= @timed sparsemixedgramcomplement(dict2, dict1, μ)
#     @show timings_stepA2[i] = t
#     col_indices = findall(reshape(nonzero_rows(ImZA),size(dict1)))
#     A1, t, _ = @timed sparseRAE(dict1, g, col_indices)
#     @show timings_stepA1[i] = t
#     A2 = ImZA[LinearIndices(size(dict1))[col_indices],:]
#     A, t, _ = @timed A1*A2
#     @show timings_stepA1A2[i] = t
#     @show nzeroels[i] = count(abs.(A.nzval) .< 1e-12)
#     @show nnzels[i] = nnz(A)
#     QRfact,t,_ = @timed qr(A)
#     @show timings_stepQR[i] = t
#     sQ[i], sT[i], sR[i], sPcol[i],sProw[i] = map(x->Base.summarysize(getfield(QRfact,x)), fieldnames(typeof(QRfact)))
# end
# @show ns
# @show timings_stepA2
# @show timings_stepA1
# @show timings_stepA1A2
# @show timings_stepQR
# @show nzeroels
# @show nnzels
# @show sQ
# @show sR
# @show sT
# @show sPcol
# @show sProw

ns = 100:100:4000
timings_stepA2 = [0.0122042, 0.0332601, 0.0632191, 0.101824, 0.150675, 0.222142, 0.268274, 0.411235, 0.426517, 0.565011, 0.66344, 0.783292, 0.977069, 1.08525, 1.17045, 1.28267, 1.47985, 1.72371, 1.95488, 2.1313, 2.26345, 2.54884, 2.98702, 3.33415, 3.29979, 3.22677, 3.49732, 3.8246, 4.07897, 4.29573, 4.56017, 4.84246, 5.11948, 5.61908, 6.06649, 6.1099, 6.65013, 7.08654, 7.30569, 7.55419]
timings_stepA1 = [0.000886431, 0.00191138, 0.00585509, 0.0047999, 0.00730389, 0.0112645, 0.0129019, 0.0197926, 0.0185817, 0.0256072, 0.0260418, 0.0332772, 0.0660462, 0.0587926, 0.0637233, 0.0660448, 0.0782921, 0.0817779, 0.0880355, 0.189295, 0.199604, 0.204784, 0.215228, 0.230622, 0.205017, 0.233195, 0.233582, 0.249021, 0.244377, 0.276348, 0.287617, 0.303176, 0.312868, 0.320966, 0.343925, 0.359172, 0.25385, 0.267193, 0.272454, 0.28739]
timings_stepA1A2 = [0.000687553, 0.00149516, 0.00245171, 0.00872468, 0.00946699, 0.0303752, 0.0167415, 0.0274834, 0.0468549, 0.148252, 0.052834, 0.17843, 0.24507, 0.271289, 0.263177, 0.262934, 0.184552, 0.194133, 0.207727, 0.227444, 0.218918, 0.261859, 0.269317, 0.338364, 0.321411, 0.356823, 0.365685, 0.407711, 0.383389, 0.410299, 0.419212, 0.441975, 0.456011, 0.574218, 0.6005, 0.613383, 0.648378, 0.671221, 0.700064, 0.71621]
timings_stepQR = [0.00559851, 0.0165652, 0.0255298, 0.122034, 0.0491162, 0.160678, 0.212403, 0.234968, 0.250166, 0.285581, 0.307682, 0.356591, 0.413038, 0.432551, 0.457001, 0.527154, 0.59226, 0.696564, 0.710568, 0.794067, 0.882694, 0.96411, 1.09587, 1.2521, 1.21601, 1.33498, 1.4257, 1.54755, 1.92415, 2.06418, 2.22292, 2.29604, 2.46212, 2.82552, 2.87405, 3.0309, 3.18702, 3.34442, 3.52284, 3.65631]
nzeroels = [32, 8, 64, 48, 72, 88, 48, 76, 84, 92, 64, 108, 0, 120, 120, 144, 74, 140, 156, 62, 184, 82, 216, 88, 92, 110, 180, 0, 94, 0, 128, 0, 0, 120, 0, 0, 0, 0, 0, 176]
nnzels = [30852, 60676, 92100, 121732, 152628, 184260, 214196, 245956, 274820, 307844, 337700, 367652, 400084, 427828, 460468, 491316, 523460, 551396, 583124, 614804, 645588, 675044, 708884, 735316, 767956, 797220, 827700, 860196, 890276, 921716, 951076, 985236, 1012788, 1045108, 1079796, 1106788, 1139524, 1166362, 1197588, 1228500]
sQ = [1083464, 2130640, 3029032, 4299960, 5143112, 6053688, 6790392, 8043288, 8678456, 9553560, 10663336, 11579624, 12733304, 13601384, 14648952, 15378584, 17354712, 17652776, 18363128, 19902432, 19696456, 22985696, 22355528, 23833128, 24782680, 25989560, 28250424, 26525656, 28006400, 29107096, 29996120, 31432184, 30984608, 33484008, 32734728, 35316392, 36597880, 37301336, 37280600, 38095832]
sR = [211848, 418856, 557736, 842808, 946712, 1110424, 1196952, 1441032, 1496504, 1669000, 1876536, 1986840, 2215688, 2381416, 2588232, 2718696, 3197864, 3199640, 3143992, 3601448, 3409464, 4271800, 3830824, 4272280, 4435560, 4442200, 5204616, 4566008, 4955928, 5129096, 5288424, 5629304, 5330968, 5958152, 5643128, 6308648, 6539704, 6562040, 6540008, 6717240]
sT = [11288, 22688, 32264, 44696, 53640, 65832, 74008, 85560, 95944, 106904, 117768, 129272, 138984, 148600, 160408, 169704, 181976, 194344, 203272, 215072, 223848, 236064, 246120, 255240, 268728, 278520, 289528, 299368, 311888, 322120, 332632, 343256, 352928, 364136, 374088, 386296, 396888, 409000, 417400, 426728]
sPcol = [5192, 10312, 15432, 20552, 25672, 30792, 35912, 41032, 46152, 51272, 56392, 61512, 66632, 71752, 76872, 81992, 87112, 92232, 97352, 102472, 107592, 112712, 117832, 122952, 128072, 133192, 138312, 143432, 148552, 153672, 158792, 163912, 169032, 174152, 179272, 184392, 189512, 194584, 199752, 204872]
sProw = [160656, 643088, 1447344, 2573280, 4021008, 5790032, 7880752, 10293872, 13027968, 16084592, 19462000, 23161488, 27182640, 31525360, 36190352, 41176848, 46485040, 52114288, 58065968, 64338768, 70933072, 77850224, 85088912, 92649104, 100530768, 108733584, 117257296, 126104464, 135272992, 144763968, 154576464, 164709424, 175164240, 185941104, 197038928, 208459056, 220202192, 232265824, 244651376, 257357936]

using PGFPlotsX, LaTeXStrings
@pgf Axis({xmode="log",ymode="log"},
    Plot({only_marks},Table(nnzels,ns.^2))
)
P1 = @pgf Axis({xlabel="nnz(A)",legend_pos="north west",legend_cell_align="left",xmode="log",ymode="log"},
    Plot({only_marks,color="blue"},Table(nnzels,timings_stepA2)),
    LegendEntry(L"I-AZ^*"),
    Plot({only_marks,color="red"},Table(nnzels,timings_stepA1)),
    LegendEntry(L"A"),
    Plot({only_marks,color="green"},Table(nnzels,timings_stepA1A2)),
    LegendEntry(L"A-AZ^*A"),
    Plot({only_marks,color="black"},Table(nnzels,timings_stepQR)),
    LegendEntry(L"QR"),
    Plot({style="black,dashed"},Table(nnzels,1e-6nnzels)),
    LegendEntry(L"\mathcal O(\textrm{nnz}(A))"),
    )
P2 = @pgf Axis({xlabel="nnz(A)",legend_pos="north west",legend_cell_align="left",xmode="log",ymode="log"},
    Plot({only_marks,color="blue"},Table(nnzels,sQ)),
    LegendEntry("factors (Q)"),
    Plot({only_marks,color="red"},Table(nnzels,sT)),
    LegendEntry(L"\tau"),
    Plot({only_marks,color="green"},Table(nnzels,sR)),
    LegendEntry(L"R"),
    Plot({only_marks,color="black"},Table(nnzels,sPcol)),
    LegendEntry("cpiv"),
    Plot({only_marks,color="brown"},Table(nnzels,sProw)),
    LegendEntry("rpivinv"),
    Plot({style="black,dashed"},Table(nnzels,nnzels)),
    LegendEntry(L"\mathcal O(\textrm{nnz}(A))"),
)



imgpath = joinpath(splitdir(@__FILE__())[1], "man", "figs")
DocumentPGFPlots.savefigs(joinpath(imgpath,"sparseAZFirstStepTimings-0"), P1)
DocumentPGFPlots.savefigs(joinpath(imgpath,"memorycheck-0"), P2)
