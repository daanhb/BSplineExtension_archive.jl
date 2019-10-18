# make number of nonzerocolumns smaller
# check correctness of difference_indices
# write nonzero_cols if domain is a grid

module AZSparse



using FrameFun.Platforms: SolverStyle
using FrameFun.BasisFunctions, FrameFun.ApproximationProblems, FrameFun.FrameFunInterface, LinearAlgebra
import FrameFun.FrameFunInterface: solver
using ..BSplineExtensionSolvers: sparseQR_solver

export AZSparseStyle
"""
    struct AZSparseStyle <: SolverStyle end

Use the sparse AZ algorithm to solve the discretized approximation problem.
"""
struct AZSparseStyle <: SolverStyle end

include("sparseAAZA.jl")

solver(style::AZSparseStyle, ap::ApproximationProblem, A::DictionaryOperator; options...) =
    solver(style, ap, A, AZ_Zt(DictionaryOperatorStyle(), ap; (options)...); options...)

function solver(::AZSparseStyle, ap::ApproximationProblem, A::AbstractOperator, Zt::AbstractOperator;
        verbose=false, REG=sparseQR_solver, threshold=default_threshold(A), options...)
    plunge_op = I-A*Zt
    dict1 = src(A)
    dict2 = dest(Zt)
    g = grid(dest(A))

    colix = nonzero_cols(dict1,g)
    verbose && @info "AZSparseStyle: create sparse A-AZ^*A"
    AAZA = sparseAAZAmatrix(dict1,dict2,g; atol=threshold)
    verbose && @info "AZSparseStyle: A-AZ^*A has size $(size(AAZA)) and $(nnz(AAZA)) nonzero elements ($(100nnz(AAZA)/prod(size(AAZA)))% fill)"
    verbose && @info "AZSparseStyle: use $(REG) as solver for first AZ step"
    psolver = IndexExtensionOperator(dict1,colix)*REG(ArrayOperator(sparseAAZAmatrix(dict1,dict2,g; atol=threshold), dict1[colix], dest(A));
        verbose=verbose, threshold=threshold, options...)
    verbose && @info "AZSparseStyle: $(REG) created"
    AZSolver(A, Zt, plunge_op, psolver)
end

end
