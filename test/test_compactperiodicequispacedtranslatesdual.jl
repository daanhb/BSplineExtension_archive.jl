

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
