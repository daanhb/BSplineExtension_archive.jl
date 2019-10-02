using StaticArrays, BSplineExtension, FrameFun, SparseArrays
using BSplineExtension.AZSparse: sparseAAZAmatrix


ns = 10:10:20
AAZA_Afactors = zeros(length(ns))
AAZA_Nfactors = zeros(length(ns))
try
    D = .4*disk() + SVector(.5,.5)
    Pbasis = NdCDBSplinePlatform((3,3))
    P = ExtensionFramePlatform(Pbasis, D)
    m = 2,2

    for (i,n) in enumerate(ns)
        N = n,n
        @show N
        dict1 = dictionary(P,N)
        g = oversampling_grid(P,N;L=m.*N)
        μ = discretemeasure(g)
        dict2 = dualdictionary(P,N,μ)

        nnzAAZA = nnz(sparseAAZAmatrix(dict1,dict2,μ))
        nnzA = nnz(sparse(AZ_A(P,N;L=m.*N)).A)
        @show AAZA_Afactors[i] = nnzAAZA/nnzA
        @show AAZA_Nfactors[i] = nnzAAZA/prod(N)
    end
finally
    @show ns
    @show AAZA_Afactors
    @show AAZA_Nfactors
end

ns = 10:5:20
AAZA_Afactors = zeros(length(ns))
AAZA_Nfactors = zeros(length(ns))
try
    D = .4*ball() + SVector(.5,.5,.5)
    Pbasis = NdCDBSplinePlatform((3,3,3))
    P = ExtensionFramePlatform(Pbasis, D)
    m = 2,2,2
    for (i,n) in enumerate(ns)
        N = n,n,n
        @show N
        dict1 = dictionary(P,N)
        g = oversampling_grid(P,N;L=m.*N)
        μ = discretemeasure(g)
        dict2 = dualdictionary(P,N,μ)

        nnzAAZA = nnz(sparseAAZAmatrix(dict1,dict2,μ))
        nnzA = nnz(sparse(AZ_A(P,N;L=m.*N)).A)
        @show AAZA_Afactors[i] = nnzAAZA/nnzA
        @show AAZA_Nfactors[i] = nnzAAZA/prod(N)
    end
finally
    @show ns
    @show AAZA_Afactors
    @show AAZA_Nfactors
end
