
using BSplineExtension, Test
using BSplineExtension: BSplineExtensionSolver
@testset "approximation power" begin
    f = (x,y)->exp(x*y)
    for d1 in 1:2, d2 in 1:2, PLATFORM in (NdEpsBSplinePlatform, NdBSplinePlatform, NdCDBSplinePlatform), N1 in (10,15), N2 in (10,15)
        P = ExtensionFramePlatform(PLATFORM((d1,d2)), (0.0..0.5)^2); N = (N1,N2)
        plunge = plungeoperator(P,N); A = AZ_A(P,N); Zt = AZ_Zt(P,N)
        M = plunge*A; S = BSplineExtensionSolver(M;crop=true)
        b = samplingoperator(P,N)*f
        x1 = S*plunge*b
        x2 = Zt*(b-A*x1)
        F = DictFun(dictionary(P,N), x1 + x2)
        @test abserror(f, F) < 1e-4
    end
end




using BSplineExtension, Test
using BSplineExtension: BSplineExtensionSolver
@testset "ExtensionFramePlatform, BSplineExtensionSolver AZ approximation power, Nd" begin
    f = (x,y)->exp(x*y)
    for d in 1:5, PLATFORM in (NdEpsBSplinePlatform, NdBSplinePlatform, NdCDBSplinePlatform)
        P = ExtensionFramePlatform(PLATFORM((d,d)), (0.0..0.5)^2); N = (30,30)
        F = Fun(f, P, N;REG=BSplineExtension.BSplineExtensionSolver, crop=true)
        @test abserror(f, F) < 1e-4
    end
end
