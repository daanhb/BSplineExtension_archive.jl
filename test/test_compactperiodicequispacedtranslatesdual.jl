

using BSplineExtension, CardinalBSplines, Test, BSplineExtension.BasisFunctions
using BSplineExtension.BasisFunctions.GridArrays: similargrid

N = 10

x = LinRange(-3,3,100)

B = GenericPeriodicEquispacedTranslates(PeriodicEquispacedGrid(N,-1.123,1.432), x->CenteredBSpline(2)(N*x),(-2..2)/N)

m = 2
D = BSplineExtension.BSplinePlatforms.CompactPeriodicEquispacedTranslatesDuals.CompactPeriodicEquispacedTranslatesDual(B, m)
m1 = evaluation_operator(B,GridArrays.similargrid(interpolation_grid(B), Float64,m*N))
m2 = evaluation_operator(D,GridArrays.similargrid(interpolation_grid(B), Float64,m*N))
@test m1'm2≈IdentityOperator(B, B)
P = platform(B)
@test P isa CDPETPlatform
μ = discretemeasure(similargrid(interpolation_grid(B),Float64,m*N))
D = dualdictionary(P, N, μ)
@test mixedgramoperator(B,D,μ)≈IdentityOperator(B, B)




P = CDPETPlatform(BSplineTranslatesBasis(6,3,-1,1))
d1 = dictionary(P,6)
d2 = azdual_dict(P,6)
g2 = mixedgramoperator(d1, d2, discretemeasure(sampling_grid(P,6)))

using InfiniteVectors
b = BSplineExtension.BSplinePlatforms.CompactPeriodicEquispacedTranslatesDuals.signal(P.dict,2)
primal_signal = PeriodicInfiniteVector(b, 12)[0:11]
c = inv(b, 2)
dual_signal = PeriodicInfiniteVector(c, 12)[0:11]

@test evaluation_operator(d1, sampling_grid(P,6)).A[:,1]≈primal_signal
@test evaluation_operator(d2, sampling_grid(P,6)).A[:,1]≈dual_signal

@test g2≈IdentityOperator(d1)
