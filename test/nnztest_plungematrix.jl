using StaticArrays, BSplineExtension, FrameFun
using BSplineExtension.AZSparse: sparseAAZAmatrix

D = .4*disk() + SVector(.5,.5)
Pbasis = NdCDBSplinePlatform((3,3))
P = ExtensionFramePlatform(Pbasis, D)
m = 2,2
ns = 10:10:1000
AAZA_Afactors = zeros(length(ns))
AAZA_Nfactors = zeros(length(ns))
for (i,n) in enumerate(ns)
    N = n,n
    dict1 = dictionary(P,N)
    g = oversampling_grid(P,N;L=m.*N)
    μ = discretemeasure(g)
    dict2 = dualdictionary(P,N,μ)

    nnzAAZA = nnz(sparseAAZAmatrix(dict1,dict2,μ))
    nnzA = nnz(sparse(AZ_A(P,N;L=m.*N)).A)
    @show AAZA_Afactors[i] = nnzAAZA/nnzA
    @show AAZA_Nfactors[i] = nnzAAZA/prod(N)
end
@show ns
@show AAZA_Afactors
@show AAZA_Nfactors
