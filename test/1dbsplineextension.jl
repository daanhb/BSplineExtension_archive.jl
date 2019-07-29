

using BSplineExtension, Test
@testset "truncated size" begin
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
    @test all(19 .<= colsizes[3,1,:,1] .<= 27)
    @test all(8 .<= colsizes[3,1,:,2] .<= 8)
    @test all(28 .<= colsizes[3,2,:,1] .<= 36)
    @test all(20 .<= colsizes[3,2,:,2] .<= 20)
    @test all(36 .<= colsizes[3,4,:,2] .<= 36)
    @test all(62 .<= colsizes[3,4,:,1][2:end] .<= 64)

    @test all(180 .<= colsizes[2,1,:,1][end-1:end] .<= 180)
    @test all(138 .<= colsizes[2,1,:,2][end-1:end] .<= 138)
    @test all(316 .<= colsizes[2,2,:,1][end-1:end] .<= 316)
    @test all(222 .<= colsizes[2,2,:,2][end-1:end] .<= 224)
    @test all(354 .<= colsizes[2,4,:,2][end-1:end] .<= 356)
    @test all(552 .<= colsizes[2,4,:,1][2:end][end-1:end] .<= 555)

    @test all(136 .<= colsizes[1,1,:,2][end-4:end] .<= 138)
    @test all(220 .<= colsizes[1,2,:,2][end-4:end] .<= 224)
    @test all(288 .<= colsizes[1,3,:,2][end-4:end] .<= 290)
    @test all(350 .<= colsizes[1,4,:,2][end-4:end] .<= 356)
end

using BSplineExtension, Test
@testset "loss of information" begin
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
@testset "errors, and timings" begin
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
