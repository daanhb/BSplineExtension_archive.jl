
using  Test
@testset "AbstractBSplinePlatforms" begin
    include("test_abstractbsplineplatforms.jl")
end

@testset "ExtensionFrame BSpline Platform" begin
    include("test_extensionframebsplineplatforms.jl")
end

@testset "Nd BSpline Platforms" begin
    include("test_ndbsplineplatforms.jl")
end

@testset "Nd ExtensionFrame BSpline Platform" begin
    include("test_ndextensionframebsplineplatforms.jl")
end

@testset "nonzero_coefficients" begin
    include("test_nonzero_coefficients.jl")
end

@testset "1D BSplineExtension" begin
    include("1dbsplineextension.jl")
end

@testset "CDPET platform" begin
    include("test_compactperiodicequispacedtranslatesdual.jl")
end

@testset "ND BSplineExtension" begin
    # include("Ndbsplineextension.jl")
end

@testset "PGFPlotsX" begin
    include("test_plots.jl")
end
