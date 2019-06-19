
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


using BSplineExtension, Test
@testset "Nd BSpline platforms, dictionaries" begin
    P = NdBSplinePlatform((1,3))
    B = dictionary(P,(10,10))

    @test FrameFun.default_discretemeasure(OversamplingStyle(),B,(11,11)) ≈ discretemeasure(PeriodicEquispacedGrid(11,0,1)^2)

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


    P = NdCDBSplinePlatform((1,3))
    d1 = dictionary(P,20)
    d2 = azdual_dict(P,20)
    g2 = mixedgramoperator(d1, d2, discretemeasure(sampling_grid(P,20)))
    @test IdentityOperator(d1)≈g2
end

using Test, BSplineExtension
@testset "samplng parameter, oversampling" begin
    P = BSplinePlatform()
    for N in 1:100
        @test divrem(samplingparameter(P, N; samplingstyle=OversamplingStyle(), oversamplingfactor=pi) ,N)[2] == 0
    end

    P = ExtensionFramePlatform(BSplinePlatform(), 0.0..0.5)
    for N in 1:100
        @test divrem(samplingparameter(P, N; samplingstyle=OversamplingStyle(), oversamplingfactor=pi) ,N)[2] == 0
    end

    P = NdBSplinePlatform((3,3))
    for N in 1:100
        @test (0,0) == map(x->divrem(x, N)[2], samplingparameter(P, (N,N); samplingstyle=ProductSamplingStyle(OversamplingStyle(),OversamplingStyle()), oversamplingfactor=pi))
    end

    P = ExtensionFramePlatform(NdBSplinePlatform((3,3)), (0.0..0.5)^2)
    for N in 5:100
        @test (0,0) == map(x->divrem(x, N)[2], samplingparameter(P, (N,N); oversamplingfactor=pi))
    end


    P = ExtensionFramePlatform(NdBSplinePlatform((1,3)),(0.0..0.5)^2)
    N = 10
    @test operator(basis(azdual_dict(P,(N,N);L=(4N,4N))))≈
        operator(basis(azdual_dict(P,(N,N),oversamplingfactor=4)))
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


using Test, BSplineExtension
@testset "BSplineTranslatesBasis, nonzero coefficients" begin
    modandsort = x->sort(mod.(x .- 1,5) .+ 1)
    for d in 1:4
        for dict in (BSplineTranslatesBasis(5,d), BSplineTranslatesBasis(10,d))
            @test modandsort(nonzero_coefficients(dict, eps())) == modandsort(findall(evaluation_matrix(dict, eps())[:] .!= 0))
            for t in .1:.1:.9
                @test modandsort(nonzero_coefficients(dict, t-eps())) == modandsort(findall(evaluation_matrix(dict, t-eps())[:] .!= 0))
                @test modandsort(nonzero_coefficients(dict, t+eps())) == modandsort(findall(evaluation_matrix(dict, t+eps())[:] .!= 0))
            end
        end
    end
end


using Test, BSplineExtension, GridArrays
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
    @test coefficient_indices_of_overlapping_elements(dict, PeriodicEquispacedGrid(100, 0,.2)) == [1,2,3,5]
    @test coefficient_indices_of_overlapping_elements(dict, PeriodicEquispacedGrid(100, .2,.4)) == [1,2,3,4]
    @test coefficient_indices_of_overlapping_elements(dict, PeriodicEquispacedGrid(100, .4,.6)) == [2,3,4,5]
    @test coefficient_indices_of_overlapping_elements(dict, PeriodicEquispacedGrid(100, .6,.8)) == [1,2,3,4,5]
    @test coefficient_indices_of_overlapping_elements(dict, PeriodicEquispacedGrid(100, .8,1)) == [1,2,4,5]
    @test coefficient_indices_of_overlapping_elements(dict, PeriodicEquispacedGrid(100, .6+eps(),.8-eps())) == [1,3,4,5]

    g = MidpointEquispacedGrid(100,UnitInterval())
    b = g[findall(GridArrays.boundary_mask(g, 0.0..0.5, true))]
    m1 = BSplineExtension.coefficient_index_mask_of_overlapping_elements(dict, b)
    m2 = sum(evaluation_matrix(dict, b);dims=1 ) .!= 0
    @test m1[:]==m2[:]
end

using BSplineExtension, Test
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

using BSplineExtension, Test
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

using BSplineExtension, Test
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



using BSplineExtension, Test
@testset "ExtensionFramePlatform, BSplineExtensionSolver truncated size" begin
    Ns = 20:20:300
    ds = 1:4
    PLATFORMs = (BSplinePlatform, EpsBSplinePlatform, CDBSplinePlatform)
    crop_tols = 10.0.^(-16.:6.:-10.)
    colsizes = Array{Int}(undef, length(PLATFORMs), length(ds), length(Ns), length(crop_tols))
    rowsizes = similar(colsizes)
    for (i,PLATFORM) in enumerate(PLATFORMs), (j,d) in enumerate(ds), (k,N) in enumerate(Ns), (l,crop_tol) in enumerate(crop_tols)
        P = ExtensionFramePlatform(PLATFORM(d), 0.0..0.5);
        plunge = plungeoperator(P,N;L=4N); A = AZ_A(P,N;L=4N); Zt = AZ_Zt(P,N;L=4N);
        M = plunge*A;
        colsizes[i,j,k,l], rowsizes[i,j,k,l]  = truncated_size(BSplineExtensionSolver(M; crop_tol=crop_tol))
    end

    @test all(rowsizes[:,1,:,:] .== 4)
    @test all(rowsizes[:,2,:,:] .== 6)
    @test all(rowsizes[:,3,:,:] .== 8)
    @test all(rowsizes[:,4,:,:] .== 10)

    # Test CDBSplinePlatform
    @test all(18 .<= colsizes[3,1,:,1] .<= 23)
    @test all(4 .<= colsizes[3,1,:,2] .<= 8)
    @test all(36 .<= colsizes[3,2,:,1] .<= 43)
    @test all(20 .<= colsizes[3,2,:,2] .<= 20)
    @test all(41 .<= colsizes[3,4,:,2] .<= 44)
    @test all(72 .<= colsizes[3,4,:,1][2:end] .<= 74)

    @test all(188 .<= colsizes[2,1,:,1][end-1:end] .<= 188)
    @test all(138 .<= colsizes[2,1,:,2][end-1:end] .<= 138)
    @test all(324 .<= colsizes[2,2,:,1][end-1:end] .<= 324)
    @test all(222 .<= colsizes[2,2,:,2][end-1:end] .<= 224)
    @test all(354 .<= colsizes[2,4,:,2][end-1:end] .<= 356)
    @test all(558 .<= colsizes[2,4,:,1][2:end][end-1:end] .<= 561)

    @test all(136 .<= colsizes[1,1,:,2][end-4:end] .<= 138)
    @test all(220 .<= colsizes[1,2,:,2][end-4:end] .<= 224)
    @test all(288 .<= colsizes[1,3,:,2][end-4:end] .<= 290)
    @test all(350 .<= colsizes[1,4,:,2][end-4:end] .<= 356)
end

using BSplineExtension, Test
@testset "ExtensionFramePlatform, BSplineExtensionSolver loss of information" begin
    Ns = [300,]
    ds = 1:4
    PLATFORMs = (BSplinePlatform, EpsBSplinePlatform, CDBSplinePlatform)
    crop_tols = [1e-8,]
    for (i,PLATFORM) in enumerate(PLATFORMs), (j,d) in enumerate(ds), (k,N) in enumerate(Ns), (l,crop_tol) in enumerate(crop_tols)
        P = ExtensionFramePlatform(PLATFORM(d), 0.0..0.5);
        plunge = plungeoperator(P,N;L=4N); A = AZ_A(P,N;L=4N); Zt = AZ_Zt(P,N;L=4N);
        M = plunge*A;
        S = BSplineExtensionSolver(M; crop_tol=crop_tol)
        @test all(size(M) .> truncated_size(S))
        @test norm(M)≈norm(S.sol.op)
    end
end


using Test, BSplineExtension, Statistics
@testset "1d spline extension approximation, errors, and timings" begin
    PLATFORMs = (EpsBSplinePlatform, BSplinePlatform, CDBSplinePlatform)
        Ns1 = [1<<k for k in 4:10]
        Ns2 = [1<<k for k in 9:16]
        ds = 1:4
        errors = Array{Float64}(undef, length(PLATFORMs), length(ds), length(Ns1))
        timings = Array{Float64}(undef, length(PLATFORMs), length(ds), length(Ns2))

    for (i,PLATFORM) in enumerate(PLATFORMs), (j,d) in enumerate(ds), (k,N) in enumerate(Ns1)
        P = ExtensionFramePlatform(PLATFORM(d), 0.0..0.5);
        F,_ = @timed Fun(exp, P, N;L=4N, REG=BSplineExtension.BSplineExtensionSolver, crop=true, crop_tol=1e-10)
        errors[i,j,k] = abserror(exp, F)
    end

    for (i,PLATFORM) in enumerate(PLATFORMs), (j,d) in enumerate(ds), (k,N) in enumerate(Ns2)
        P = ExtensionFramePlatform(PLATFORM(d), 0.0..0.5);
        timings[i,j,k]= median([@timed(Fun(exp, P, N;L=4N, REG=BSplineExtension.BSplineExtensionSolver, crop=true, crop_tol=1e-10))[2] for l in 1:4])
    end

    # Test if errors go down fast enough
    for (i,PLATFORM) in enumerate(PLATFORMs), (j,d) in enumerate(ds)
        @test all(errors[i,j,:] .<= max.(Float64.(Ns1).^(-d)*errors[i,j,1]*Ns1[1]^d, 1e-10))
    end

    # Test if method is fast enough
    for (i,PLATFORM) in enumerate(PLATFORMs), (j,d) in enumerate(ds)
        @test all(timings[i,j,:] .<= 5e-5Ns2)
    end

end
