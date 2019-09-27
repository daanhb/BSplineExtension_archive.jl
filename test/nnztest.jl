using Pkg;Pkg.activate(localprojectdir())
using BSplineExtension, AMD, Metis, SparseArrays, LinearAlgebra, StaticArrays

Pbasis = NdCDBSplinePlatform((3,3,3))
ns = [16,32,40,45,50]
nnz1 = zeros(Int,length(ns))
nnz2 = similar(nnz1)
nnz3 = similar(nnz1)

D = .4ball() + SVector(.5,.5,.5)
P = ExtensionFramePlatform(Pbasis, D)



N = (16,16,16)
@show N
Abasis = AZ_A(Pbasis,N;oversamplingfactor=2^3,samplingstyle=ProductSamplingStyle([OversamplingStyle() for i in 1:3]...))
A = AZ_A(P,N;samplingstyle=OversamplingStyle())
Zt = AZ_Zt(P,N;samplingstyle=OversamplingStyle())

M = A-A*Zt*A
m = Matrix(M)
s = sparse(M).A
s.colptr
s.nzval
s.rowval
nzcols = BSplineExtension.BSplineExtensionSolvers.nonzero_cols(basis(src(M)), supergrid(grid(dest(M))), D)
length(basis(src(M)))
count(abs.(s.nzval) .<= 1e-12)
s_red = s[:,LinearIndices(size(basis(src(M))))[nzcols]]
count(abs.(s_red.nzval) .<= 1e-12)
s_res = copy(s); s_res[:,LinearIndices(size(basis(src(M))))[nzcols]] .= 0; dropzeros!(s_res)
nnz(s_res)/prod(size(s))
maximum(abs.(s_res))


for (i,n) in enumerate(ns)
    N = ntuple(k->n,Val(3))
    @show N
    Abasis = AZ_A(Pbasis,N;oversamplingfactor=2^3,samplingstyle=ProductSamplingStyle([OversamplingStyle() for i in 1:3]...))
    A = AZ_A(P,N;oversamplingfactor=2^3,samplingstyle=ProductSamplingStyle([OversamplingStyle() for i in 1:3]...))
    Zt = AZ_A(P,N;oversamplingfactor=2^3,samplingstyle=ProductSamplingStyle([OversamplingStyle() for i in 1:3]...))
    Asparse = sparse(Abasis'Abasis).A
    nnz1[i] = nnz(Asparse)
    @show nnz1[i]
    @time perm, iperm = Metis.permutation(Asparse)
    @time nnz2[i] = nnz(cholesky(Asparse[perm,perm]) )
    @show nnz2[i]
    @time p = amd(Asparse)
    @time nnz3[i] = nnz(cholesky(Asparse[p,p]) )
    @show nnz3[i]
end
@show nnz1, nnz2, nnz3


# stevin:test vincentcp$ julia nnztest.jl
# N = (16, 16, 16)
# nnz1[i] = 1404928
#   0.179578 seconds (156.28 k allocations: 13.231 MiB)
#   1.275051 seconds (1.80 M allocations: 429.267 MiB, 17.17% gc time)
# nnz2[i] = 5701160
#   0.229743 seconds (419.40 k allocations: 44.969 MiB, 34.39% gc time)
#   0.483440 seconds (506.98 k allocations: 327.677 MiB, 3.29% gc time)
# nnz3[i] = 5160345
# N = (32, 32, 32)
# nnz1[i] = 11239424
#   1.429207 seconds (18 allocations: 43.501 MiB, 8.78% gc time)
#  29.582258 seconds (90 allocations: 10.118 GiB, 0.98% gc time)
# nnz2[i] = 250985473
#   0.595304 seconds (16 allocations: 191.101 MiB, 1.28% gc time)
#  18.270003 seconds (86 allocations: 7.015 GiB, 0.55% gc time)
# nnz3[i] = 189172914
# N = (40, 40, 40)
# nnz1[i] = 46656000
#   4.582240 seconds (18 allocations: 179.200 MiB, 4.39% gc time)
# 307.849967 seconds (89 allocations: 41.576 GiB, 0.36% gc time)
# nnz2[i] = 1045077502
#   4.195013 seconds (16 allocations: 787.891 MiB, 1.88% gc time)
# 188.150797 seconds (86 allocations: 32.961 GiB, 0.19% gc time)
# nnz3[i] = 839160690
# N = (45, 45, 45)
# nnz1[i] = 66430125
#   7.113253 seconds (18 allocations: 255.150 MiB, 3.14% gc time)
# 1009.818471 seconds (89 allocations: 76.203 GiB, 0.14% gc time)
# nnz2[i] = 1935561063
#   6.410364 seconds (16 allocations: 1.096 GiB, 9.53% gc time)
# 977.556854 seconds (86 allocations: 71.386 GiB, 0.09% gc time)
# nnz3[i] = 1785527518
