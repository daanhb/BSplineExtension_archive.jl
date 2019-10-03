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
        REG=sparseQR_solver, threshold=default_threshold(A), options...)
    plunge_op = I-A*Zt
    dict1 = src(A)
    dict2 = dest(Zt)
    g = grid(dest(A))

    colix = nonzero_cols(dict1,g)

    AAZA = sparseAAZAmatrix(dict1,dict2,g; atol=threshold)
    psolver = IndexExtensionOperator(dict1,colix)*REG(ArrayOperator(sparseAAZAmatrix(dict1,dict2,g; atol=threshold), dict1[colix], dest(A)); threshold=threshold, options...)
    AZSolver(A, Zt, plunge_op, psolver)
end

end
