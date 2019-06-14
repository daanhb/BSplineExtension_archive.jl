
using LinearAlgebra, BSplineExtension, Test
@testset "BSpline platforms, dual dictionaries" begin


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
end
