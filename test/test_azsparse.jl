using StaticArrays, BSplineExtension, FrameFun, LinearAlgebra, SparseArrays, BasisFunctions, Test


D = .4*disk() + SVector(.5,.5)
Pbasis = NdCDBSplinePlatform((1,1))
P1 = ExtensionFramePlatform(Pbasis, D)
N1 = 20,20
m1 = (2,2)

D = .4*disk() + SVector(.5,.5)
Pbasis = NdCDBSplinePlatform((3,3))
P3 = ExtensionFramePlatform(Pbasis, D)
N3 = 20,20
m3 = (2,2)

D = .4*disk() + SVector(.5,.5)
Pbasis = NdCDBSplinePlatform((2,2))
P2 = ExtensionFramePlatform(Pbasis, D)
N2 = 20,20
m2 = (2,3)

D = .4*ball() + SVector(.5,.5,.5)
Pbasis = NdCDBSplinePlatform((2,2,2))
P4 = ExtensionFramePlatform(Pbasis, D)
N4 = 20,20,10
m4 = (2,2,2)

using BSplineExtension.AZSparse: bsplinesignal, difference_indices, sparsemixedgramcomplement_nzband, sparsemixedgramcomplement
@testset "AZSparse: sparsemixedgramcomplement_nzband, sparsemixedgramcomplement" begin
    for (P,N,m) in zip((P1,P2,P3,P4), (N1,N2,N3,N4), (m1,m2,m3,m4))
        dict1 = dictionary(P,N)
        g = oversampling_grid(P,N;L=m.*N)
        γ = supergrid(g)
        μ = discretemeasure(g)
        dict2 = dualdictionary(P,N,μ)
        L = LinearIndices(size(dict1))
        C = CartesianIndices(size(dict1))

        G = Matrix(mixedgramoperator(dict2,dict1,μ))
        indices = nonzero_cols(dict1, μ)
        GG = copy(G);GG[:,L[indices]] .= 0;
        @test sum(abs.(GG))≈count(sum(abs.(GG),dims=1) .≈1 )

        IG = I-G
        IGE = IG[:,L[indices]]

        b = bsplinesignal(dict1,μ)
        b̃ = bsplinesignal(dict2,μ)

        dff = difference_indices(b,b̃,m)

        R = Matrix{Float64}(undef, length(dff), length(indices))

        M = GridArrays.ModCartesianIndicesBase.ModCartesianIndices(size(basis(dict1)),first(dff), last(dff))
        M1 = GridArrays.ModCartesianIndicesBase.ModCartesianIndices(size(basis(dict1)),CartesianIndex(ntuple(k->1,Val(dimension(dict1)))),CartesianIndex(size(basis(dict1))))

        for (i,k) in enumerate(indices)
            kdff = k .+ dff
            M = GridArrays.ModCartesianIndicesBase.ModCartesianIndices(size(basis(dict1)),first(kdff),last(kdff))
            for (j,l) in enumerate(M)
                R[j,i] = IG[L[M1[l]],L[k]]
            end
        end
        @test sum(abs.(R),dims=1) ≈ sum(abs.(IGE),dims=1)
        @test count(sum(abs.(R),dims=2) .+ 1 .≈ 1) == 0

        RR = sparsemixedgramcomplement_nzband(indices, Float64, b̃, support(b̃), b, support(b), m, dff, g)
        @test RR≈R

        @test nnz(sparse(RR))<= nnz(sparse(R))
        s = sparsemixedgramcomplement(dict2,dict1,μ)
        @test IGE≈s
    end
end
