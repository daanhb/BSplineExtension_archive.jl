
using PGFPlotsX
import PGFPlotsX: Plot, PlotInc, Options, Axis

for F in (:DictFun, :Expansion)
    @eval begin
        Axis(F::$F, trailing...; opts...) =
            Axis(Options(), F, trailing...; opts...)
        function Axis(options::Options, F::$F, trailing...; opts...)
            if !(haskey(options, :ymode))
                @pgf  options = {options..., ymode="log"}
            end
            Axis(options, F, trailing...; opts...)
        end
    end
end

for plot in (:Plot, :PlotInc)
    _plot = Meta.parse("_"*string(plot))
    @eval begin
        $(plot)(F::FrameFun.DictFun; opts...) =
            $(plot)(Options(), F; opts...)

        function $(plot)(options::Options, F::DictFun; plot_extension=false, opts...)
            $(plot)(options, plot_extension ? Expansion(basis(F),coefficients(F)) : expansion(F); opts...)
        end

        $(plot)(f::Function, F::DictFun; opts...) =
            $(plot)(Options(), f, F; opts...)
        $(plot)(F::DictFun, f::Function; opts...) =
            $(plot)(Options(), f, F)
        $(plot)(F::DictFun, f::DictFun; opts...) =
            $(plot)(Options(), f, F; opts...)

        $(plot)(options::Options, F::DictFun, f::Function; opts...) =
            $(_plot)(options, f, F; opts...)
        $(plot)(options::Options, f::Function, F::DictFun; opts...) =
            $(_plot)(options, f, F; opts...)
        $(plot)(options::Options, F::DictFun, f::DictFun; opts...) =
            $(_plot)(options, f, F; opts...)
        function $(_plot)(options::Options, f::Function, F::DictFun; plot_extension=false, opts...)
            if plot_extension
                $(plot)(options, Expansion(basis(F),coefficients(F)), f; opts...)
            else
                $(plot)(options, expansion(F), f; opts...)
            end
        end

        function $(plot)(options::Options, F::Expansion; n=200, plot_complex=false, opts...)
            grid = plotgrid(dictionary(F), n)
            vals = plot_complex ? F(grid) : real.(F(grid))
            options = @pgf {options..., no_markers}
            $(plot)(options, Table([grid, BasisFunctions.postprocess(dictionary(F), grid, vals)]))
        end

        $(plot)(options::Options, F::Expansion, f::Function; opts...) =
            $(_plot)(options, f, F; opts...)

        $(plot)(options::Options, f::Function, F::Expansion; opts...) =
            $(_plot)(options, f, F; opts...)

        function $(_plot)(options::Options, f::Function, F::Expansion; plot_complex=false, n=200, opts...)
            grid = plotgrid(dictionary(F), n)
            vals = abs.(f.(grid) - F.(grid))
            options = @pgf {options..., no_markers}
            $(plot)(options, Table([grid, BasisFunctions.postprocess(dictionary(F), grid, vals)]))
        end
    end
end
