

using Test, BSplineExtension
@testset "BSplineTranslatesBasis, nonzero coefficients" begin
    modandsort = x->sort(mod.(x .- 1,5) .+ 1)
    for d in 1:4
        for dict in (BSplineTranslatesBasis(5,d), BSplineTranslatesBasis(10,d))
            @test modandsort(nonzero_coefficients(dict, eps())[1]) == modandsort(findall(evaluation_matrix(dict, eps())[:] .!= 0))
            for t in .1:.1:.9
                @test modandsort(nonzero_coefficients(dict, t-eps())[1]) == modandsort(findall(evaluation_matrix(dict, t-eps())[:] .!= 0))
                @test modandsort(nonzero_coefficients(dict, t+eps())[1]) == modandsort(findall(evaluation_matrix(dict, t+eps())[:] .!= 0))
            end
        end
    end

    modandsort2 = x->(sort(mod.(x[1] .- 1,5) .+ 1), sort(mod.(x[2] .- 1,5) .+ 1))
    function getnonzeroindices(dict, x)
        m = abs.(reshape(evaluation_matrix(dict, [x]),size(dict)))
        modandsort(findall(sum(m;dims=2)[:] .!=0)), modandsort(findall(sum(m;dims=1)[:] .!=0))
    end

    for d1 in 1:4, d2 in 1:4
        dict = dictionary(NdBSplinePlatform((d1,d2)),(5,5))
        modandsort2(nonzero_coefficients(dict, (eps(),eps()))) == getnonzeroindices(dict, (eps(),eps()))
        for t1 in .1:.1:.9, t2 in .1:.1:.9
            for x in ((t1-eps(),t2-eps()),(t1+eps(),t2-eps()),(t1-eps(),t2+eps()),(t1+eps(),t2+eps()))
                @test modandsort2(nonzero_coefficients(dict, x)) == getnonzeroindices(dict, x)
            end
        end
    end
end


using Test, BSplineExtension, GridArrays
using BSplineExtension.BSplineExtensionSolvers: coefficient_indices_of_overlapping_elements,
    coefficient_index_mask_of_overlapping_elements
@testset "Extensionframe platform, nonzero indices" begin
    dict = BSplineTranslatesBasis(5,2)
    @test coefficient_indices_of_overlapping_elements(dict, PeriodicEquispacedGrid(100, 0,.1)) == [1,2,5]
    @test coefficient_indices_of_overlapping_elements(dict, PeriodicEquispacedGrid(100, .1,.3)) == [1,2,3]
    @test coefficient_indices_of_overlapping_elements(dict, PeriodicEquispacedGrid(100, .3,.5)) == [1,2,3,4]
    @test coefficient_indices_of_overlapping_elements(dict, PeriodicEquispacedGrid(100, .5,.7)) == [3,4,5]
    @test coefficient_indices_of_overlapping_elements(dict, PeriodicEquispacedGrid(100, .7,.9)) == [1,3,4,5]
    @test coefficient_indices_of_overlapping_elements(dict, PeriodicEquispacedGrid(100, .9,1.)) == [1,2,5]

    @test coefficient_indices_of_overlapping_elements(dict, PeriodicEquispacedGrid(100, 0+eps(),.1-eps())) == [1,2,5]
    @test coefficient_indices_of_overlapping_elements(dict, PeriodicEquispacedGrid(100, .1+eps(),.3-eps())) == [1,2,3]
    @test coefficient_indices_of_overlapping_elements(dict, PeriodicEquispacedGrid(100, .3+eps(),.5-eps())) == [2,3,4]
    @test coefficient_indices_of_overlapping_elements(dict, PeriodicEquispacedGrid(100, .5+eps(),.7-eps())) == [3,4,5]
    @test coefficient_indices_of_overlapping_elements(dict, PeriodicEquispacedGrid(100, .7+eps(),.9-eps())) == [1,4,5]
    @test coefficient_indices_of_overlapping_elements(dict, PeriodicEquispacedGrid(100, .9+eps(),1. -eps())) == [1,2,5]

    dict = BSplineTranslatesBasis(5,3)
    @test coefficient_indices_of_overlapping_elements(dict, PeriodicEquispacedGrid(100, 0,.2)) == [1,2,3,5]
    @test coefficient_indices_of_overlapping_elements(dict, PeriodicEquispacedGrid(100, .2,.4)) == [1,2,3,4]
    @test coefficient_indices_of_overlapping_elements(dict, PeriodicEquispacedGrid(100, .4,.6)) == [2,3,4,5]
    @test coefficient_indices_of_overlapping_elements(dict, PeriodicEquispacedGrid(100, .6,.8)) == [1,2,3,4,5]
    @test coefficient_indices_of_overlapping_elements(dict, PeriodicEquispacedGrid(100, .8,1)) == [1,2,4,5]
    @test coefficient_indices_of_overlapping_elements(dict, PeriodicEquispacedGrid(100, .6+eps(),.8-eps())) == [1,3,4,5]

    g = MidpointEquispacedGrid(100,UnitInterval())
    b = g[findall(GridArrays.boundary_mask(g, 0.0..0.5, true))]
    m1 = coefficient_index_mask_of_overlapping_elements(dict, b)
    m2 = sum(evaluation_matrix(dict, b);dims=1 ) .!= 0
    @test m1[:]==m2[:]
end

using BSplineExtension, Test
using BSplineExtension.BSplineExtensionSolvers: nonzero_rows, nonzero_cols
@testset "ExtensionFramePlatform, nonzero_rows, nonzero_cols" begin
    # 1D nonzero_cols
    for d in (1,2,3,4)
        P = ExtensionFramePlatform(EpsBSplinePlatform(d), 0.0..0.5)
        N = 100
        A = Matrix(AZ_A(P,N; L=4N))
        Zt = Matrix(AZ_Zt(P,N;L=4N))
        gt = nonzero_cols(basis(dictionary(P,N)), supergrid(sampling_grid(P,N;L=4N)), 0.0..0.5)
        n = findall(sum(abs.(Matrix(A*Zt*A-A)); dims=1)[:] .> 1e-10)
        for ni in n
            @test ni ∈ gt
        end
    end

    # 1D nonzero_rows
    P = FrameFun.ExtensionFramePlatform(CDBSplinePlatform(), 0.0..0.5)
    N = 100
    M = plungeoperator(P,N;L=4N)*AZ_A(P,N;L=4N);size(M)
    @test findall(nonzero_rows(Matrix(M),nonzero_tol=1e-10))[:] ==
        findall(sum(abs.(Matrix(M));dims=2)[:] .> 1e-10)

    # 2D nonzero_cols
    P = ExtensionFramePlatform(NdEpsBSplinePlatform((3,3)),(0.0..0.5)^2)
    N = (10,10)
    A = Matrix(AZ_A(P,N))
    Zt = Matrix(AZ_Zt(P,N))
    gt = nonzero_cols(basis(dictionary(P,N)), supergrid(sampling_grid(P,N)), (0.0..0.5)^2)
    n = findall(reshape(sum(abs.(Matrix(A*Zt*A-A)); dims=1),N) .> 1e-10)
    for ni in n
        @test ni ∈ gt
    end

    # 2D nonzero_rows
    P = ExtensionFramePlatform(NdCDBSplinePlatform((3,3)),(0.0..0.5)^2)
    N = (10,10)
    M = plungeoperator(P,N)*AZ_A(P,N);size(M)
    @test findall(nonzero_rows(Matrix(M),nonzero_tol=1e-10))[:] ==
        findall(sum(abs.(Matrix(M));dims=2)[:] .> 1e-10)
end
