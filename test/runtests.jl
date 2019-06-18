
using LinearAlgebra, BSplineExtension, Test
using BSplineExtension.FrameFun: approximationproblem

@testset "BSpline platforms, dictionaries" begin


    P = BSplinePlatform()
    B = dictionary(P,10)

    @test FrameFun.default_discretemeasure(OversamplingStyle(),B,11) ≈ discretemeasure(PeriodicEquispacedGrid(11,0,1))
    ap = FrameFun.approximationproblem(P,10)
    μ1 = FrameFun.discretemeasure(SamplingStyle(ap), ap)
    μ2 = discretemeasure(sampling_grid(P,10))
    @test μ1≈μ2


    d1 = dictionary(P,10)
    d2 = azdual_dict(P,10)

    g1 = mixedgramoperator(d1, d2)
    g2 = mixedgramoperator(d1, d2, discretemeasure(sampling_grid(P,10)))
    @test g1 isa CirculantOperator
    @test g2 isa CirculantOperator
    @test g2 ≈ IdentityOperator(B)

    P = EpsBSplinePlatform()
    g = sampling_grid(P,10)
    d1 = dictionary(P,1000)
    d2 = azdual_dict(P,1000;threshold=1e-4)
    g2 = mixedgramoperator(d1, d2, discretemeasure(sampling_grid(P,1000)))
    @test norm(IdentityOperator(d1)-g2) < 1e-3


    P = CDBSplinePlatform(5)
    d1 = dictionary(P,20)
    d2 = azdual_dict(P,20)
    g2 = mixedgramoperator(d1, d2, discretemeasure(sampling_grid(P,20)))
    @test IdentityOperator(d1)≈g2


    P = EpsBSplinePlatform()
    @test 20==length(sampling_grid(P,10; oversamplingfactor=1.6,samplingstyle=OversamplingStyle()))
    @test 20==FrameFun.samplingparameter(P,10; oversamplingfactor=1.6,samplingstyle=OversamplingStyle())

end

using Test, BSplineExtension
@testset "Basis platform approximation" begin

    P1 = BSplinePlatform()
    P2 = EpsBSplinePlatform()
    P3 = CDBSplinePlatform()
    @test SamplingStyle(approximationproblem(P1,10)) == InterpolationStyle()
    @test SamplingStyle(approximationproblem(P2,10)) == InterpolationStyle()
    @test SamplingStyle(approximationproblem(P3,10)) == OversamplingStyle()

    # Test the approximation power
    f = x->exp(cos(10pi*x))
    N = 60
    x = PeriodicEquispacedGrid(10N, support(dictionary(P1,10)))
    t = .123
    for P in (P1,P2,P3)
        @test SolverStyle(InterpolationStyle(), approximationproblem(P1, 10)) == DualStyle()
        F = Fun(f, P, N)
        @test norm(F.(x)-f.(x),Inf)  < .005
        @test abs(F(t)-f(t)) < .0002
    end
end

using Test, BSplineExtension, GridArrays, DomainSets
using BSplineExtension: coefficient_indices_of_overlapping_elements

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
    @test coefficient_indices_of_overlapping_elements(dict, PeriodicEquispacedGrid(100, 0,.2)) == [1,2,4,5]
    @test coefficient_indices_of_overlapping_elements(dict, PeriodicEquispacedGrid(100, .2,.4)) == [1,2,3,5]
    @test coefficient_indices_of_overlapping_elements(dict, PeriodicEquispacedGrid(100, .4,.6)) == [1,2,3,4]
    @test coefficient_indices_of_overlapping_elements(dict, PeriodicEquispacedGrid(100, .6,.8)) == [1,2,3,4,5]
    @test coefficient_indices_of_overlapping_elements(dict, PeriodicEquispacedGrid(100, .8,1)) == [1,3,4,5]
    @test coefficient_indices_of_overlapping_elements(dict, PeriodicEquispacedGrid(100, .6+eps(),.8-eps())) == [2,3,4,5]

    g = MidpointEquispacedGrid(100,UnitInterval())
    b = g[findall(GridArrays.boundary_mask(g, 0.0..0.5, true))]
    m1 = BSplineExtension.coefficient_index_mask_of_overlapping_elements(dict, b)
    m2 = sum(evaluation_matrix(dict, b);dims=1 ) .!= 0
    @test m1[:]==m2[:]
end

using BSplineExtension, Test, DomainSets
using BSplineExtension.FrameFun: ExtensionFramePlatform
using BSplineExtension: nonzero_rows, nonzero_cols
@testset "ExtensionFramePlatform, nonzero_rows, nonzero_cols" begin
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

    P = FrameFun.ExtensionFramePlatform(CDBSplinePlatform(), 0.0..0.5)
    N = 100
    M = plungeoperator(P,N;L=4N)*AZ_A(P,N;L=4N);size(M)
    @test findall(nonzero_rows(Matrix(M),nonzero_tol=1e-10))[:] ==
        findall(sum(abs.(Matrix(M));dims=2)[:] .> 1e-10)
end

using BSplineExtension, Test, DomainSets
using BSplineExtension.FrameFun: ExtensionFramePlatform
using BSplineExtension: BSplineExtensionSolver
@testset "ExtensionFramePlatform, BSplineExtensionSolver approximation power" begin
        for d in 1:5
            for PLATFORM in (EpsBSplinePlatform, BSplinePlatform, CDBSplinePlatform)
                P = ExtensionFramePlatform(PLATFORM(d), 0.0..0.5); N = 30
                plunge = plungeoperator(P,N;L=4N); A = AZ_A(P,N;L=4N); Zt = AZ_Zt(P,N;L=4N)
                M = plunge*A; S = BSplineExtensionSolver(M;crop=true)
                b = samplingoperator(P,N;L=4N)*exp
                x1 = S*plunge*b
                x2 = Zt*(b-A*x1)
                F = DictFun(dictionary(P,N), x1 + x2)
                @test abserror(exp, F) < 1e-4
            end
        end
end

using BSplineExtension, Test, DomainSets
using BSplineExtension.FrameFun: ExtensionFramePlatform
using BSplineExtension: BSplineExtensionSolver
@testset "ExtensionFramePlatform, BSplineExtensionSolver AZ approximation power" begin
        for d in 1:5
            for PLATFORM in (EpsBSplinePlatform, BSplinePlatform, CDBSplinePlatform)
                P = ExtensionFramePlatform(EpsBSplinePlatform(1), 0.0..0.5); N = 30
                F = Fun(exp, P, N;REG=BSplineExtension.BSplineExtensionSolver, L=4N, crop=true)
                @test abserror(exp, F) < 1e-4
            end
        end
end
