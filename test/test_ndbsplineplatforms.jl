

using BSplineExtension, Test
@testset "(dual)dictionaries" begin
    P = NdBSplinePlatform((1,3))
    B = dictionary(P,(10,10))


    d1 = dictionary(P,(10,10))
    d2 = azdual_dict(P,(10,10))

    g1 = mixedgramoperator(d1, d2)
    sampling_grid(P,(10,10))
    g2 = mixedgramoperator(d1, d2, discretemeasure(sampling_grid(P,(10,10))))
    all(isa.(elements(g1), CirculantOperator))
    @test all(isa.(elements(g1), CirculantOperator))
    @test all(isa.(elements(g2), CirculantOperator))
    @test g2 ≈ IdentityOperator(B)

    P = NdEpsBSplinePlatform((1,3))
    g = sampling_grid(P,10)
    d1 = dictionary(P,30)
    d2 = azdual_dict(P,30;threshold=1e-4)
    g2 = mixedgramoperator(d1, d2, discretemeasure(sampling_grid(P,30)))
    @test norm(IdentityOperator(d1)-g2) < 1e-3

    opts = (oversamplingfactor=4,)
    P = NdCDBSplinePlatform((1,3))
    d1 = dictionary(P,20)
    d2 = azdual_dict(P,20; opts...)
    g2 = mixedgramoperator(d1, d2, discretemeasure(sampling_grid(P,20; opts...)))
    @test IdentityOperator(d1)≈g2
end


using BSplineExtension, Test
@testset "approximation" begin
    opts = (oversamplingfactor=4,)
    P1 = NdBSplinePlatform((3,3))
    P2 = NdEpsBSplinePlatform((3,3))
    P3 = NdCDBSplinePlatform((3,3))
    f = (x,y) -> exp(cos(10pi*sqrt((x-.5)^2+(y-.5)^2)))

    N = (100,100)
    x = PeriodicEquispacedGrid(3N[1], UnitInterval())^2
    t = (.123,.203)
    for P in (P1,P2,P3)
        @test SolverStyle(ProductSamplingStyle(InterpolationStyle(),InterpolationStyle()), approximationproblem(P, (10,10))) ==
            ProductSolverStyle(DualStyle(),DualStyle())
        F = Fun(f, P, N; (elements(SamplingStyle(P))[1]==InterpolationStyle() ? tuple() : opts )...)
        @test norm(F(x)-[f(xi...) for xi in x],Inf)  < .1
        @test abs(F(t...)-f(t...)) < 1e-5
    end
end
