using StaticArrays, BSplineExtension, FrameFun, LinearAlgebra, SparseArrays, BasisFunctions, Test


D = .4*disk() + SVector(.5,.5)
    Pbasis = NdCDBSplinePlatform((1,1))
    P1 = ExtensionFramePlatform(Pbasis, D)
    N1 = 20,20
    m1 = (2,2)
    f1 = (x,y)->exp(x*y)

    D = .4*disk() + SVector(.5,.5)
    Pbasis = NdCDBSplinePlatform((3,3))
    P3 = ExtensionFramePlatform(Pbasis, D)
    N3 = 20,20
    m3 = (2,2)
    f3 = (x,y)->exp(x*y)

    D = .4*disk() + SVector(.5,.5)
    Pbasis = NdCDBSplinePlatform((2,2))
    P2 = ExtensionFramePlatform(Pbasis, D)
    N2 = 20,20
    m2 = (2,3)
    f2 = (x,y)->exp(x*y)

    D = .4*ball() + SVector(.5,.5,.5)
    Pbasis = NdCDBSplinePlatform((2,2,2))
    P4 = ExtensionFramePlatform(Pbasis, D)
    N4 = 20,20,10
    m4 = (2,2,2)
    f4 = (x,y,z)->exp(x*y*z)

    D = .4*ball() + SVector(.5,.5,.5)
    Pbasis = NdCDBSplinePlatform((3,3,3))
    P5 = ExtensionFramePlatform(Pbasis, D)
    N5 = 20,20,30
    m5 = (2,2,2)
    f5 = (x,y,z)->exp(x*y*z)

using BSplineExtension.AZSparse: compactinfinitevector, difference_indices, sparsemixedgramcomplement_nzband, sparsemixedgramcomplement
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

        b = compactinfinitevector(dict1,g)
        b̃ = compactinfinitevector(dict2,g)

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

using BSplineExtension.AZSparse: sparsemixedgramcomplement, sparseRAE
using BSplineExtension.BSplineExtensionSolvers: nonzero_rows
using SparseArrays, Test
@testset "AZSparse: sparseRAE" begin
    for (P,N,m) in zip((P1,P2,P3,P4), (N1,N2,N3,N4), (m1,m2,m3,m4))
        dict1 = dictionary(P,N)
        g = oversampling_grid(P,N;L=m.*N)
        γ = supergrid(g)
        μ = discretemeasure(g)
        dict2 = dualdictionary(P,N,μ)

        s = sparsemixedgramcomplement(dict2,dict1,μ)
        col_indices = findall(reshape(nonzero_rows(s),size(dict1)))
        RAE =sparseRAE(dict1,g,col_indices)
        RAEref = sparse((E = IndexExtensionOperator(dict1, col_indices);tmp=AZ_A(P,N;L=m.*N)*E;tmp=Matrix(tmp);tmp[abs.(tmp).<1e-14].=0;tmp))
        @test RAE≈RAEref
    end
end

using BSplineExtension.AZSparse: sparseAAZAmatrix, sparseidentity
@testset "AZSparse: sparseAAZAmatrix" begin
    for (P,N,m) in zip((P1,P2,P3,P4), (N1,N2,N3,N4), (m1,m2,m3,m4))
        dict1 = dictionary(P,N)
        g = oversampling_grid(P,N;L=m.*N)
        γ = supergrid(g)
        μ = discretemeasure(g)
        dict2 = dualdictionary(P,N,μ)
        L = LinearIndices(size(dict1))

        AAZA = sparseAAZAmatrix(dict1,dict2,g)
        ImZA = sparsemixedgramcomplement(dict2,dict1,μ)

        A_ref = Matrix(AZ_A(P,N;L=m.*N))
        Z_ref = Matrix(AZ_Z(P,N;L=m.*N))
        ix1 = nonzero_cols(dict1, g)
        A = sparseRAE(dict1, g, ix1)
        ix2 = nonzero_cols(dict2, g)
        Z = sparseRAE(dict2, g, ix2)

        @test Z'A ≈ (Z_ref'A_ref)[L[ix2],L[ix1]]
        I_ref = sparse(Matrix{Int}(I, length(dict1),length(dict1))[L[ix2],L[ix1]] )

        II = sparseidentity(ix1,ix2)
        @test I_ref ≈ II
        @test I_ref - Z'A ≈ (I-Z_ref'A_ref)[L[ix2],L[ix1]]
        indices = LinearIndices(size(dict1))[nonzero_cols(dict1,μ)]
        @test (LinearAlgebra.I-Z_ref'*A_ref)[:,indices] ≈ ImZA
        AAZAref = A_ref-A_ref*Z_ref'*A_ref
        @test (AAZAref)[:,indices] ≈ AAZA
        @test (AAZAref)[:,indices]≈AAZA
        (AAZAref)[:,indices].=0
        @test norm(AAZAref,1) <1e-9
    end
end

@testset "AZSparse: FrameFunInterface" begin
    err = [1e-2,1e-2,1e-6,1e-1]
    for (P,N,m,f,e) in zip((P1,P2,P3,P4), (N1,N2,N3,N4), (m1,m2,m3,m4), (f1,f2,f3,f4), err)
        F = Fun(f, P, N; L=m.*N,solverstyle=AZSparseStyle())
        @test abserror(f,F) < e
    end
end
