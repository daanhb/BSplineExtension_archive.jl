using BSplineExtension.AZSparse: sparseRAE, sparseidentity
using DomainSets, BSplineExtension, StaticArrays, LinearAlgebra, SparseArrays
D = .4*disk() + SVector(.5,.5)
Pbasis = NdCDBSplinePlatform((2,2))
P = ExtensionFramePlatform(Pbasis, D)
N = 20,20
m = (2,2)
ns = 500:500:6000

n = 100
N = (n,n)
dict1 = dictionary(P,N)
g = oversampling_grid(P,N;L=m.*N)
μ = discretemeasure(g)
dict2 = dualdictionary(P,N,μ)

ix1, t, _ = @timed nonzero_cols(dict1, g)
A, t, _ = @timed sparseRAE(dict1, g, ix1)
ix2, t, _ = @timed nonzero_cols(dict2, g)
Z, t, _ = @timed sparseRAE(dict2, g, ix2)
ImZA, t, _ = @timed sparseidentity(ix1,ix2)-Z'A
RAE, t, _ = @timed sparseRAE(dict1, g, ix2)
M, t, _ = @timed RAE*ImZA
AA, t, _ = @timed sparse(AZ_A(P,N;L=m.*N)).A
# QRfactA,t,_ = @timed qr(AA);
QRfact,t,_ = @timed qr(M);
Base.summarysize(QRfact.Q) + Base.summarysize(QRfact.R) + Base.summarysize(QRfact.pcol) + Base.summarysize(QRfact.prow)
Base.summarysize(QRfact)

timings_stepA = zeros(length(ns))
timings_stepix1 = similar(timings_stepA)
timings_stepix2 = similar(timings_stepA)
timings_stepZ = similar(timings_stepA)
timings_stepImZA = similar(timings_stepA)
timings_stepRAE = similar(timings_stepA)
timings_stepM = similar(timings_stepA)
timings_stepAA = similar(timings_stepA)
timings_stepQRM = similar(timings_stepA)
timings_stepQRA = similar(timings_stepA)
nzeroels = Vector{Int}(undef,length(ns))
nnzelsM = Vector{Int}(undef,length(ns))
nnzelsA = Vector{Int}(undef,length(ns))
sQ = Vector{Int}(undef,length(ns))
sR = Vector{Int}(undef,length(ns))
sT = Vector{Int}(undef,length(ns))
sProw = Vector{Int}(undef,length(ns))
sPcol = Vector{Int}(undef,length(ns))
for (i,n) in enumerate(ns)
    N = (n,n)
    @show N
    dict1 = dictionary(P,N)
    g = oversampling_grid(P,N;L=m.*N)
    μ = discretemeasure(g)
    dict2 = dualdictionary(P,N,μ)


    ix1, t, _ = @timed nonzero_cols(dict1, g)
    @show timings_stepix1[i] = t
    A, t, _ = @timed sparseRAE(dict1, g, ix1)
    @show timings_stepA[i] = t
    ix2, t, _ = @timed nonzero_cols(dict2, g)
    @show timings_stepix2[i] = t
    Z, t, _ = @timed sparseRAE(dict2, g, ix2)
    @show timings_stepZ[i] = t
    ImZA, t, _ = @timed sparseidentity(ix1,ix2)-Z'A
    @show timings_stepImZA[i] = t
    RAE, t, _ = @timed sparseRAE(dict1, g, ix2)
    @show timings_stepRAE[i] = t
    M, t, _ = @timed RAE*ImZA
    @show timings_stepM[i] = t

    AA, t, _ = @timed sparse(AZ_A(P,N;L=m.*N)).A
    @show timings_stepAA[i] = t

    @show nzeroels[i] = count(abs.(M.nzval) .< 1e-12)
    @show nnzelsM[i] = nnz(M)
    @show nnzelsA[i] = nnz(AA)

    QRfactM,t,_ = @timed qr(M);
    @show timings_stepQRM[i] = t
    # QRfactA,t,_ = @timed qr(AA);
    # @show timings_stepQRA[i] = t
    sQ[i], sT[i], sR[i], sPcol[i],sProw[i] = map(x->Base.summarysize(getfield(QRfactM,x)), fieldnames(typeof(QRfactM)))
end
@show ns
@show timings_stepix1
@show timings_stepA
@show timings_stepix2
@show timings_stepZ
@show timings_stepImZA
@show timings_stepRAE
@show timings_stepM
@show timings_stepQRA
@show timings_stepQRM
@show timings_stepAA
@show nzeroels
@show nnzelsA
@show nnzelsM
@show sQ
@show sR
@show sT
@show sPcol
@show sProw

# ns = 100:100:4000
# timings_stepA2 = [0.0122042, 0.0332601, 0.0632191, 0.101824, 0.150675, 0.222142, 0.268274, 0.411235, 0.426517, 0.565011, 0.66344, 0.783292, 0.977069, 1.08525, 1.17045, 1.28267, 1.47985, 1.72371, 1.95488, 2.1313, 2.26345, 2.54884, 2.98702, 3.33415, 3.29979, 3.22677, 3.49732, 3.8246, 4.07897, 4.29573, 4.56017, 4.84246, 5.11948, 5.61908, 6.06649, 6.1099, 6.65013, 7.08654, 7.30569, 7.55419]
# timings_stepA1 = [0.000886431, 0.00191138, 0.00585509, 0.0047999, 0.00730389, 0.0112645, 0.0129019, 0.0197926, 0.0185817, 0.0256072, 0.0260418, 0.0332772, 0.0660462, 0.0587926, 0.0637233, 0.0660448, 0.0782921, 0.0817779, 0.0880355, 0.189295, 0.199604, 0.204784, 0.215228, 0.230622, 0.205017, 0.233195, 0.233582, 0.249021, 0.244377, 0.276348, 0.287617, 0.303176, 0.312868, 0.320966, 0.343925, 0.359172, 0.25385, 0.267193, 0.272454, 0.28739]
# timings_stepA1A2 = [0.000687553, 0.00149516, 0.00245171, 0.00872468, 0.00946699, 0.0303752, 0.0167415, 0.0274834, 0.0468549, 0.148252, 0.052834, 0.17843, 0.24507, 0.271289, 0.263177, 0.262934, 0.184552, 0.194133, 0.207727, 0.227444, 0.218918, 0.261859, 0.269317, 0.338364, 0.321411, 0.356823, 0.365685, 0.407711, 0.383389, 0.410299, 0.419212, 0.441975, 0.456011, 0.574218, 0.6005, 0.613383, 0.648378, 0.671221, 0.700064, 0.71621]
# timings_stepQR = [0.00559851, 0.0165652, 0.0255298, 0.122034, 0.0491162, 0.160678, 0.212403, 0.234968, 0.250166, 0.285581, 0.307682, 0.356591, 0.413038, 0.432551, 0.457001, 0.527154, 0.59226, 0.696564, 0.710568, 0.794067, 0.882694, 0.96411, 1.09587, 1.2521, 1.21601, 1.33498, 1.4257, 1.54755, 1.92415, 2.06418, 2.22292, 2.29604, 2.46212, 2.82552, 2.87405, 3.0309, 3.18702, 3.34442, 3.52284, 3.65631]
# nzeroels = [32, 8, 64, 48, 72, 88, 48, 76, 84, 92, 64, 108, 0, 120, 120, 144, 74, 140, 156, 62, 184, 82, 216, 88, 92, 110, 180, 0, 94, 0, 128, 0, 0, 120, 0, 0, 0, 0, 0, 176]
# nnzels = [30852, 60676, 92100, 121732, 152628, 184260, 214196, 245956, 274820, 307844, 337700, 367652, 400084, 427828, 460468, 491316, 523460, 551396, 583124, 614804, 645588, 675044, 708884, 735316, 767956, 797220, 827700, 860196, 890276, 921716, 951076, 985236, 1012788, 1045108, 1079796, 1106788, 1139524, 1166362, 1197588, 1228500]
# sQ = [1083464, 2130640, 3029032, 4299960, 5143112, 6053688, 6790392, 8043288, 8678456, 9553560, 10663336, 11579624, 12733304, 13601384, 14648952, 15378584, 17354712, 17652776, 18363128, 19902432, 19696456, 22985696, 22355528, 23833128, 24782680, 25989560, 28250424, 26525656, 28006400, 29107096, 29996120, 31432184, 30984608, 33484008, 32734728, 35316392, 36597880, 37301336, 37280600, 38095832]
# sR = [211848, 418856, 557736, 842808, 946712, 1110424, 1196952, 1441032, 1496504, 1669000, 1876536, 1986840, 2215688, 2381416, 2588232, 2718696, 3197864, 3199640, 3143992, 3601448, 3409464, 4271800, 3830824, 4272280, 4435560, 4442200, 5204616, 4566008, 4955928, 5129096, 5288424, 5629304, 5330968, 5958152, 5643128, 6308648, 6539704, 6562040, 6540008, 6717240]
# sT = [11288, 22688, 32264, 44696, 53640, 65832, 74008, 85560, 95944, 106904, 117768, 129272, 138984, 148600, 160408, 169704, 181976, 194344, 203272, 215072, 223848, 236064, 246120, 255240, 268728, 278520, 289528, 299368, 311888, 322120, 332632, 343256, 352928, 364136, 374088, 386296, 396888, 409000, 417400, 426728]
# sPcol = [5192, 10312, 15432, 20552, 25672, 30792, 35912, 41032, 46152, 51272, 56392, 61512, 66632, 71752, 76872, 81992, 87112, 92232, 97352, 102472, 107592, 112712, 117832, 122952, 128072, 133192, 138312, 143432, 148552, 153672, 158792, 163912, 169032, 174152, 179272, 184392, 189512, 194584, 199752, 204872]
# sProw = [160656, 643088, 1447344, 2573280, 4021008, 5790032, 7880752, 10293872, 13027968, 16084592, 19462000, 23161488, 27182640, 31525360, 36190352, 41176848, 46485040, 52114288, 58065968, 64338768, 70933072, 77850224, 85088912, 92649104, 100530768, 108733584, 117257296, 126104464, 135272992, 144763968, 154576464, 164709424, 175164240, 185941104, 197038928, 208459056, 220202192, 232265824, 244651376, 257357936]
using PGFPlotsX, LaTeXStrings
# @pgf Axis({xmode="log",ymode="log"},
#     Plot({only_marks},Table(nnzels,ns.^2))
# )
# P1 = @pgf Axis({xlabel="nnz(A)",legend_pos="north west",legend_cell_align="left",xmode="log",ymode="log"},
#     Plot({only_marks,color="blue"},Table(nnzels,timings_stepA2)),
#     LegendEntry(L"I-AZ^*"),
#     Plot({only_marks,color="red"},Table(nnzels,timings_stepA1)),
#     LegendEntry(L"A"),
#     Plot({only_marks,color="green"},Table(nnzels,timings_stepA1A2)),
#     LegendEntry(L"A-AZ^*A"),
#     Plot({only_marks,color="black"},Table(nnzels,timings_stepQR)),
#     LegendEntry(L"QR"),
#     Plot({style="black,dashed"},Table(nnzels,1e-6nnzels)),
#     LegendEntry(L"\mathcal O(\textrm{nnz}(A))"),
#     )
# P2 = @pgf Axis({xlabel="nnz(A)",legend_pos="north west",legend_cell_align="left",xmode="log",ymode="log"},
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
#     LegendEntry(L"\mathcal O(\textrm{nnz}(A))"),
# )
#
#
#
# imgpath = joinpath(splitdir(@__FILE__())[1], "man", "figs")
# DocumentPGFPlots.savefigs(joinpath(imgpath,"sparseAZFirstStepTimings-0"), P1)
# DocumentPGFPlots.savefigs(joinpath(imgpath,"memorycheck-0"), P2)
#
#
#
# ns = 100:100:5000
# timings_stepA2 = [0.0176167, 0.0561136, 0.0859398, 0.146345, 0.205741, 0.281015, 0.371867, 0.478616, 0.592617, 0.724643, 0.846592, 0.993884, 1.20849, 1.36286, 1.55354, 1.72371, 1.92928, 1.97149, 2.17903, 2.39969, 2.61558, 2.86883, 3.11978, 3.39086, 3.67131, 4.01034, 4.73229, 5.12963, 5.46009, 5.8163, 6.20691, 6.5958, 7.20956, 7.62395, 7.72699, 7.72372, 8.07545, 8.48133, 8.96708, 10.4021, 11.3011, 11.7851, 12.4593, 12.7531, 13.0539, 13.4306, 14.025, 14.585, 14.5761, 14.3472]
# timings_stepA1 = [0.00107623, 0.00269999, 0.00468025, 0.0069396, 0.00945233, 0.0108746, 0.0184074, 0.0178873, 0.0196143, 0.0190919, 0.0383438, 0.0401187, 0.0793121, 0.0741334, 0.0941964, 0.0988057, 0.106399, 0.112502, 0.188632, 0.196552, 0.205603, 0.223209, 0.23339, 0.233118, 0.225374, 0.233648, 0.314951, 0.356248, 0.329451, 0.398273, 0.36846, 0.410787, 0.463215, 0.433678, 0.35985, 0.277666, 0.284846, 0.30168, 0.313358, 0.405134, 0.454005, 0.491446, 0.472941, 0.451444, 0.528745, 0.54727, 0.531313, 0.641644, 0.552451, 0.563873]
# timings_stepA1A2 = [0.000849201, 0.00192841, 0.00325465, 0.00929494, 0.00927682, 0.0139946, 0.0151184, 0.0188761, 0.0602746, 0.175656, 0.117749, 0.239781, 0.337376, 0.335149, 0.325926, 0.417567, 0.252639, 0.213026, 0.238058, 0.248325, 0.25969, 0.277467, 0.288466, 0.300168, 0.31611, 0.348157, 0.460426, 0.474729, 0.482681, 0.454627, 0.477316, 0.521126, 0.617785, 0.716005, 0.549754, 0.580317, 0.585552, 0.615107, 0.631675, 0.923239, 0.93902, 1.05659, 0.999335, 0.882925, 0.766151, 0.806422, 0.826149, 1.00071, 0.814309, 0.822821]
# timings_stepQR = [0.00710678, 0.0731015, 0.231087, 0.0563421, 0.251521, 0.278437, 0.303349, 0.342536, 0.339492, 0.470254, 0.538709, 0.591916, 0.547393, 0.569933, 0.642947, 0.975678, 0.82541, 0.887265, 0.960487, 1.06305, 1.14574, 1.2509, 1.35394, 1.47232, 1.54907, 2.18743, 2.57096, 2.74135, 2.80853, 3.08647, 3.0947, 3.43551, 3.66343, 3.23757, 2.87407, 3.0885, 3.21865, 3.35438, 3.57454, 5.72183, 5.56133, 5.88952, 6.25944, 4.85481, 5.00003, 5.15482, 5.4269, 6.40318, 5.53977, 5.72688]
# nzeroels = [32, 8, 64, 48, 72, 88, 48, 76, 112, 92, 64, 108, 0, 120, 120, 72, 0, 140, 156, 62, 92, 82, 216, 0, 92, 110, 180, 0, 188, 0, 128, 0, 0, 0, 0, 0, 0, 0, 0, 0, 680, 0, 0, 0, 0, 0, 0, 0, 0, 208]
# nnzels = [30852, 60676, 92100, 121732, 152628, 184260, 214196, 245956, 274820, 307844, 337700, 367652, 400084, 427828, 460468, 491316, 523460, 551396, 583124, 614804, 645588, 675044, 708884, 735316, 767956, 797220, 827700, 860196, 890276, 921716, 951076, 985236, 1012788, 1045108, 1079796, 1106788, 1139524, 1166362, 1197588, 1228500, 1260836, 1286932, 1323524, 1352644, 1379556, 1416244, 1442948, 1479700, 1502714, 1541188]
# sQ = [1083048, 2128496, 3112840, 4311064, 5144360, 6034728, 6799880, 7797112, 8499800, 9315624, 10898648, 11584072, 12807192, 14375288, 14484168, 15988568, 17344248, 18157064, 18334632, 19865680, 20244696, 23034704, 22644264, 22506568, 23592872, 25873880, 26401016, 27181080, 27390224, 29119928, 29191304, 32659736, 31715696, 33427016, 33906072, 35376792, 36429464, 36561016, 36502632, 37312680, 40143656, 40313176, 41963720, 41469528, 42640152, 44929448, 45191688, 45453096, 46193272, 48912720]
# sR = [211816, 418888, 582392, 845480, 947512, 1108808, 1169368, 1372296, 1478312, 1624280, 1969688, 2034680, 2239944, 2613400, 2567544, 2891544, 3198072, 3279608, 3205944, 3590072, 3575992, 4281720, 4015464, 3881768, 4074680, 4603240, 4642056, 4766520, 4766664, 5137432, 5048248, 5958824, 5547608, 5938824, 5991336, 6322856, 6484056, 6449912, 6313592, 6477960, 7115848, 6967464, 7372664, 7177048, 7421896, 7875608, 7897880, 7861176, 8012184, 8656360]
# sT = [11288, 22688, 32264, 44696, 53640, 65832, 74008, 85560, 95944, 106904, 117768, 129272, 138984, 148600, 160408, 169704, 181976, 194344, 203272, 215072, 223848, 236064, 246120, 255240, 268728, 278520, 289528, 299368, 311888, 322120, 332632, 343256, 352928, 364136, 374088, 386296, 396888, 409000, 417400, 426728, 440264, 451240, 459992, 472296, 483816, 492728, 504296, 514424, 525976, 537024]
# sPcol = [5192, 10312, 15432, 20552, 25672, 30792, 35912, 41032, 46152, 51272, 56392, 61512, 66632, 71752, 76872, 81992, 87112, 92232, 97352, 102472, 107592, 112712, 117832, 122952, 128072, 133192, 138312, 143432, 148552, 153672, 158792, 163912, 169032, 174152, 179272, 184392, 189512, 194584, 199752, 204872, 209992, 215112, 220232, 225352, 230472, 235592, 240712, 245832, 250904, 256072]
# sProw = [160656, 643088, 1447344, 2573280, 4021008, 5790032, 7880752, 10293872, 13027968, 16084592, 19462000, 23161488, 27182640, 31525360, 36190352, 41176848, 46485040, 52114288, 58065968, 64338768, 70933072, 77850224, 85088912, 92649104, 100530768, 108733584, 117257296, 126104464, 135272992, 144763968, 154576464, 164709424, 175164240, 185941104, 197038928, 208459056, 220202192, 232265824, 244651376, 257357936, 270387344, 283736880, 297409648, 311402768, 325718576, 340356800, 355316144, 370596016, 386198640, 402122640]
# P1 = @pgf Axis({xlabel="nnz(A)",legend_pos="north west",legend_cell_align="left",xmode="log",ymode="log"},
#     Plot({only_marks,color="blue"},Table(nnzels,timings_stepA2)),
#     LegendEntry(L"I-AZ^*"),
#     Plot({only_marks,color="red"},Table(nnzels,timings_stepA1)),
#     LegendEntry(L"A"),
#     Plot({only_marks,color="green"},Table(nnzels,timings_stepA1A2)),
#     LegendEntry(L"A-AZ^*A"),
#     Plot({only_marks,color="black"},Table(nnzels,timings_stepQR)),
#     LegendEntry(L"QR"),
#     Plot({style="black,dashed"},Table(nnzels,1e-6nnzels)),
#     LegendEntry(L"\mathcal O(\textrm{nnz}(A))"),
#     )
# P2 = @pgf Axis({xlabel="nnz(A)",legend_pos="north west",legend_cell_align="left",xmode="log",ymode="log"},
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
#     LegendEntry(L"\mathcal O(\textrm{nnz}(A))"),
# )


# ns = 100:100:4000
# timings_stepix1 = [0.004469089, 0.017844914, 0.040683492, 0.070782912, 0.11107062, 0.160514051, 0.218331554, 0.291812075, 0.359433935, 0.444151013, 0.536974026, 0.638629377, 0.751264882, 0.874124103, 1.03523213, 1.151381234, 1.289985313, 1.443624802, 1.606352014, 1.78694054, 1.973633376, 2.162271281, 2.362101437, 2.560390724, 2.789595531, 3.01669268, 3.253071005, 3.499921165, 3.752968459, 4.279813684, 4.385442173, 4.611978552, 4.860191108, 5.32392887, 5.482752119, 6.026568829, 6.373535962, 6.667447352, 7.269612105, 7.458315073]
# timings_stepA = [0.000651179, 0.001643786, 0.002717255, 0.005219634, 0.006254944, 0.009090422, 0.009813687, 0.013636778, 0.016441298, 0.019184116, 0.023134706, 0.027024272, 0.031917117, 0.035262829, 0.040058695, 0.045365372, 0.051037896, 0.054528539, 0.060648198, 0.068795872, 0.079569076, 0.083564563, 0.092067956, 0.099869899, 0.105870999, 0.117638157, 0.128747162, 0.147426381, 0.16072609, 0.153260869, 0.165887634, 0.172709831, 0.186644544, 0.210485351, 0.225891599, 0.248567269, 0.301228778, 0.270744924, 0.27544376, 0.299305362]
# timings_stepix2 = [0.008225859, 0.03455366, 0.076385685, 0.13709728, 0.214969253, 0.310126643, 0.420972879, 0.551971629, 0.702368796, 0.864089644, 1.051221428, 1.242019219, 1.469250065, 1.741540045, 1.970246985, 2.285867425, 2.508488562, 2.810302644, 3.144521693, 3.48475979, 3.83019535, 4.183590522, 4.590698957, 4.998892045, 5.412202962, 6.083068421, 6.328680859, 6.812484869, 7.314860976, 7.849557593, 8.353607581, 8.918477612, 9.467999338, 10.233716346, 10.69439121, 11.999361401, 12.707192791, 13.080729336, 13.649870322, 14.518970289]
# timings_stepZ = [0.001537527, 0.003249425, 0.004991674, 0.007494897, 0.008022552, 0.012098554, 0.013137195, 0.037266206, 0.024847826, 0.030441218, 0.03564074, 0.040575526, 0.062592045, 0.091307104, 0.07181858, 0.093925885, 0.091429343, 0.093194268, 0.101446264, 0.107581077, 0.115685528, 0.121440565, 0.129650966, 0.141955174, 0.249170496, 0.309433027, 0.268515743, 0.278649251, 0.311955369, 0.414549025, 0.328137525, 0.347762662, 0.357411616, 0.370150194, 0.389768352, 0.525911198, 0.44631148, 0.447877673, 0.477449507, 0.478094788]
# timings_stepImZA = [0.000511412, 0.001393634, 0.095005149, 0.002554275, 0.004062554, 0.007899495, 0.007448651, 0.009212162, 0.010307385, 0.032786148, 0.031833749, 0.033472546, 0.116639302, 0.160133871, 0.125733406, 0.144734712, 0.137839303, 0.142756226, 0.144258323, 0.16463561, 0.144464521, 0.143680652, 0.153734769, 0.162750572, 0.065556343, 0.098687392, 0.084148174, 0.080427945, 0.1002211, 0.198613568, 0.112126002, 0.11854677, 0.124366549, 0.220808087, 0.232894288, 0.339352069, 0.242733769, 0.241800574, 0.266479469, 0.272389806]
# timings_stepRAE = [0.000914977, 0.001997877, 0.003333009, 0.006426486, 0.008233602, 0.005812338, 0.023500984, 0.103242089, 0.032882231, 0.111658666, 0.113391976, 0.118487031, 0.03838817, 0.043748867, 0.046329104, 0.052947505, 0.056801057, 0.063771713, 0.069186999, 0.075429644, 0.083925204, 0.088139479, 0.095228753, 0.106485562, 0.121739902, 0.153634514, 0.132740561, 0.144801139, 0.140189691, 0.241702062, 0.161727902, 0.173169711, 0.180722149, 0.191654214, 0.203058868, 0.271818266, 0.247628979, 0.249277348, 0.270195794, 0.287427434]
# timings_stepM = [0.000936195, 0.002341148, 0.003221057, 0.002942604, 0.031197355, 0.012931682, 0.112915947, 0.021989079, 0.122008958, 0.035104091, 0.03973522, 0.061230783, 0.065606432, 0.107996653, 0.088326748, 0.089030457, 0.098590828, 0.109794356, 0.121501114, 0.125281748, 0.13357657, 0.164562846, 0.172487873, 0.185188624, 0.201184562, 0.270571747, 0.237652648, 0.242668875, 0.414074919, 0.583653837, 0.455775019, 0.471702883, 0.501596808, 0.52080691, 0.530976556, 0.618717076, 0.59408438, 0.772952927, 0.755067766, 0.816627105]
# timings_stepQR = [0.003214203, 0.028105823, 0.019000856, 0.025439552, 0.131817458, 0.170342207, 0.072225391, 0.20108118, 0.222356155, 0.254856397, 0.274793236, 0.297218122, 0.345795661, 0.518991127, 0.431525078, 0.484904878, 0.533096096, 0.580935114, 0.631147619, 0.871367492, 0.782346894, 0.931578658, 1.003241485, 1.094173653, 1.205095669, 1.410533778, 1.507980268, 1.599239517, 1.999779465, 2.347308494, 2.142324097, 2.255603804, 2.417723245, 2.504436429, 3.040889582, 3.596390084, 2.960271104, 3.326774087, 3.574361091, 3.466489267]
# nzeroels = [11254, 26030, 35372, 49396, 62351, 76824, 87304, 90082, 111626, 92715, 129884, 127282, 133594, 156880, 162840, 148878, 149887, 183566, 209328, 201066, 203724, 207978, 247180, 217392, 235583, 220348, 196950, 245448, 240271, 290368, 260246, 275436, 302244, 290242, 313024, 301016, 282850, 324322, 296742, 310498]
# nnzels = [7852, 16188, 24028, 32588, 40524, 48300, 56380, 64428, 72908, 80348, 87756, 96988, 103916, 113500, 120604, 128844, 136060, 145388, 152748, 160236, 168652, 177356, 183868, 193980, 200604, 209980, 217308, 226028, 232876, 241628, 249852, 257196, 265772, 272588, 279644, 289004, 296332, 305178, 313916, 321612]
# sQ = [465672, 886408, 1432968, 1835768, 2282360, 2759520, 3212816, 3756552, 4151760, 4421968, 4799968, 5443312, 5903944, 6574840, 6583992, 7121616, 7502184, 7994736, 8558400, 8909656, 9399248, 10005840, 10251648, 10631216, 11057680, 11654904, 11910440, 12568112, 13196128, 13372024, 13924232, 14120968, 14787720, 15177408, 15504736, 16079832, 16403656, 16932552, 17291384, 17803816]
# sR = [101976, 192488, 295736, 406680, 502504, 600952, 697432, 800600, 889144, 933112, 1021928, 1156904, 1283432, 1386568, 1404584, 1516808, 1592552, 1712776, 1829688, 1935176, 2023496, 2150376, 2222696, 2230792, 2327064, 2508760, 2517704, 2627176, 2835960, 2850008, 2991672, 3019256, 3117320, 3207576, 3318984, 3430536, 3488872, 3570712, 3694376, 3810648]
# sT = [7144, 14232, 21832, 30504, 37352, 45200, 52720, 59688, 68320, 73632, 82176, 91104, 98232, 105912, 109672, 118608, 127544, 132384, 139200, 148424, 154720, 163920, 171360, 176448, 182464, 190984, 196776, 203984, 217088, 220568, 226792, 231576, 244312, 250400, 257280, 264312, 271624, 281064, 285976, 294680]
# sPcol = [5192, 10312, 15432, 20552, 25672, 30792, 35912, 41032, 46152, 51272, 56392, 61512, 66632, 71752, 76872, 81992, 87112, 92232, 97352, 102472, 107592, 112712, 117832, 122952, 128072, 133192, 138312, 143432, 148552, 153672, 158792, 163912, 169032, 174152, 179272, 184392, 189512, 194584, 199752, 204872]
# sProw = [160656, 643088, 1447344, 2573280, 4021008, 5790032, 7880752, 10293872, 13027968, 16084592, 19462000, 23161488, 27182640, 31525360, 36190352, 41176848, 46485040, 52114288, 58065968, 64338768, 70933072, 77850224, 85088912, 92649104, 100530768, 108733584, 117257296, 126104464, 135272992, 144763968, 154576464, 164709424, 175164240, 185941104, 197038928, 208459056, 220202192, 232265824, 244651376, 257357936]



ns = 500:500:6000
timings_stepix1 = [0.11087458, 0.459462904, 1.032279905, 1.79900089, 2.828245927, 4.048425776, 5.526741985, 7.210118497, 9.167677836, 11.431338812, 13.908526684, 16.652517453]
timings_stepA = [0.005853802, 0.018859633, 0.040772106, 0.087102194, 0.128501669, 0.174912339, 0.330866938, 0.409790278, 0.462656331, 0.53563116, 0.619866546, 0.737325228]
timings_stepix2 = [0.214860609, 0.882706181, 1.992691909, 3.507229297, 5.629427266, 7.891436061, 10.837623655, 14.07484646, 17.710945415, 21.972876289, 26.788492296, 32.366217676]
timings_stepZ = [0.014130249, 0.026583246, 0.056893186, 0.214501328, 0.263599473, 0.322673618, 0.264141763, 0.334387721, 0.536405298, 0.611890028, 0.661083822, 0.785226583]
timings_stepImZA = [0.003506995, 0.014334141, 0.04334511, 0.049242071, 0.146638985, 0.207285726, 0.257737783, 0.295971817, 0.254872633, 0.284821024, 0.340188947, 0.567436922]
timings_stepRAE = [0.007120244, 0.022689097, 0.13959216, 0.188746886, 0.120776883, 0.154523898, 0.213887351, 0.275286776, 0.387082754, 0.455566338, 0.636715779, 0.572025572]
timings_stepM = [0.014128477, 0.038347616, 0.095738689, 0.128355292, 0.21623179, 0.322649332, 0.528696383, 0.791080101, 0.905108352, 1.046714309, 1.133489346, 1.412069959]
timings_stepQRA = [2.3202927533e-314, 2.316879139e-314, 2.320292793e-314, 2.320292714e-314, 5.0e-324, 2.3202929114e-314, 2.320292951e-314, 2.32029303e-314, 2.57e-322, 2.37e-322, 3.0e-323, 2.2380424363e-314]
timings_stepQRM = [0.207681471, 0.436487419, 0.78856212, 1.219242721, 1.900954518, 2.836989422, 4.663170607, 4.909185936, 7.776115759, 9.171055877, 9.308484157, 11.340587151]
timings_stepAA = [0.559428242, 2.367553163, 4.282126464, 6.752661379, 10.523343063, 18.456038405, 36.429539369, 51.627128706, 77.024307569, 115.188544896, 126.092840764, 172.931701586]
nzeroels = [62351, 92715, 162840, 201066, 235583, 290368, 313024, 310498, 338160, 461648, 320289, 478694]
nnzelsA = [4523589, 24629573, 55416805, 72381069, 113097069, 162859419, 301715757, 289527633, 498757589, 615750477, 547388901, 886681513]
nnzelsM = [214907, 427107, 661378, 877164, 1086886, 1325602, 1551092, 1729106, 1962550, 2214680, 2354153, 2590632]
sQ = [7596792, 15281824, 22974888, 30057352, 35848464, 44065984, 49135736, 54373192, 57914008, 66598296, 70405248, 74619856]
sR = [1584392, 3197144, 4810744, 6150184, 7223432, 8925944, 9697336, 10730632, 10878936, 12593608, 13168952, 13686456]
sT = [73896, 148096, 222632, 293368, 366784, 446512, 520440, 601608, 660536, 747192, 807760, 878544]
sPcol = [25672, 51272, 76872, 102472, 128072, 153672, 179272, 204872, 230472, 256072, 281672, 307272]
sProw = [4021008, 16084592, 36190352, 64338768, 100530768, 144763968, 197038928, 257357936, 325718576, 402122640, 486567952, 579056976]

P1 = @pgf Axis({xlabel="nnz(A)",legend_pos="north west",legend_cell_align="left",xmode="log",ymode="log"},
    PlotInc({},Table(nnzelsM,timings_stepix1)),
    LegendEntry(L"ix1"),
    PlotInc({},Table(nnzelsM,timings_stepA)),
    LegendEntry(L"A"),
    PlotInc({},Table(nnzelsM,timings_stepix2)),
    LegendEntry(L"ix2"),
    PlotInc({},Table(nnzelsM,timings_stepZ)),
    LegendEntry(L"Z"),
    PlotInc({},Table(nnzelsM,timings_stepImZA)),
    LegendEntry(L"ImZA"),
    PlotInc({},Table(nnzelsM,timings_stepRAE)),
    LegendEntry(L"RAE"),
    PlotInc({},Table(nnzelsM,timings_stepM)),
    LegendEntry(L"M"),
    PlotInc({},Table(nnzelsM,timings_stepAA)),
    LegendEntry(L"A"),
    # PlotInc({},Table(nnzelsM,timings_stepQRA)),
    # LegendEntry(L"QR_A"),
    PlotInc({},Table(nnzelsM,timings_stepQRM)),
    LegendEntry(L"QR_M"),
    )


P2 = @pgf Axis({xlabel="nnz(A)",legend_pos="north west",legend_cell_align="left",xmode="log",ymode="log"},
    Plot({only_marks,color="blue"},Table(nnzelsM,sQ)),
    LegendEntry("factors (Q)"),
    Plot({only_marks,color="red"},Table(nnzelsM,sT)),
    LegendEntry(L"\tau"),
    Plot({only_marks,color="green"},Table(nnzelsM,sR)),
    LegendEntry(L"R"),
    Plot({only_marks,color="black"},Table(nnzelsM,sPcol)),
    LegendEntry("cpiv"),
    Plot({only_marks,color="brown"},Table(nnzelsM,sProw)),
    LegendEntry("rpivinv"),
    Plot({style="black,dashed"},Table(nnzelsM,nnzelsM)),
    LegendEntry(L"\mathcal O(\textrm{nnz}(A))"),
)

P3 = @pgf Axis({xlabel="nnz(A)",legend_pos="north west",legend_cell_align="left",xmode="log",ymode="log"},
    Plot({only_marks,color="blue"},Table(nnzelsM,nnzelsA)),
    LegendEntry("factors (Q)"),
    Plot({only_marks,color="red"},Table(nnzelsM,nnzelsM)),
)
