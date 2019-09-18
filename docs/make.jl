using Pkg
Pkg.develop(PackageSpec(path=splitdir(@__DIR__)[1]))
using Documenter, BSplineExtension, BSplineExtension.BSplinePlatforms, BSplineExtension.BSplineExtensionSolvers


include("render_figs.jl")


const render_pdf = "pdf" in ARGS
let r = r"buildroot=(.+)", i = findfirst(x -> occursin(r, x), ARGS)
    global const buildroot = i === nothing ? (@__DIR__) : first(match(r, ARGS[i]).captures)
end

const format = if render_pdf
    LaTeX(
        platform = "texplatform=docker" in ARGS ? "docker" : "native"
    )
else
    Documenter.HTML(
        prettyurls = ("deploy" in ARGS),
    )
end
ENV["JULIA_DEBUG"] = ""
makedocs(sitename="BSplineExtension.jl",
    modules = [BSplineExtension,BSplinePlatforms,BSplineExtensionSolvers],
    authors = "vincentcp",
    format = format,
    pages = [
        "Home" => "index.md",
        "Manual" => ["Basis platform" => "man/basisplatform.md",
                    "1D extension platform" => "man/1dextension.md",
                    "ND extension platform" => "man/ndbsplineextension.md"]
        ],
    doctest=true
)

if "deploy" in ARGS && Sys.ARCH === :x86_64 && Sys.KERNEL === :Linux
    deploydocs(
        repo = "github.com/FrameFunVC/BSplineExtension.jl.git",
    )
end
