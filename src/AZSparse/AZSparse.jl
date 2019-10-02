# make number of nonzerocolumns smaller
# check correctness of difference_indices
# write nonzero_cols if domain is a grid

module AZSparse



using FrameFun.Platforms: SolverStyle
using FrameFun.BasisFunctions, FrameFun.ApproximationProblems, FrameFun.FrameFunInterface
import FrameFun.FrameFunInterface: solver
"""
    struct AZSparseStyle <: SolverStyle end

Use the sparse AZ algorithm to solve the discretized approximation problem.
"""
struct AZSparseStyle <: SolverStyle end

include("sparseAAZA.jl")
# struct SparseAZSolver{T} <: VectorizingSolverOperator{T}
#     op      :: DictionaryOperator{T}
#     grid_res     :: DictionaryOperator{T}
#     dict_ext     :: DictionaryOperator{T}
#     sol     :: VectorizingSolverOperator{T}
#
#     dict_scratch
#     grid_scratch
#
#     src_linear  ::  Vector{T}
#     dest_linear ::  Vector{T}
#
# end
#
# function linearized_apply!(op::SparseAZSolver, dest::Vector, src::Vector)
#     apply!(op.grid_res, op.grid_scratch, src)
#     apply!(op.sol, op.dict_scratch, op.grid_scratch)
#     apply!(op.dict_ext, dest, op.dict_scratch)
# end
#
# apply!(s::SparseAZSolver, coef_dest, coef_src) = _apply!(s, coef_dest, coef_src,
#         s.plunge_op, s.A, s.Zt, s.psolver, s.Pb, s.x1)
#
# function _apply!(s::SparseAZSolver, coef_dest, coef_src, plunge_op::DictionaryOperator, A, Zt, psolver, Pb, x1)
#     # Step 1: Solve (I-A*Zt)*A * x1 = (I-A*Zt)*b
#     apply!(plunge_op, Pb, coef_src)
#     apply!(psolver, x1, Pb)
#
#     # Step 2: Compute x2 =  Zt*(b-A*x1)
#     apply!(A, Pb, x1)
#     Pb .= coef_src .- Pb
#     apply!(Zt, coef_dest, Pb)
#
#     # Step 3: x = x1 + x2
#     coef_dest .+= x1
# end

function solver(::AZSparseStyle, ap::ApproximationProblem, A::AbstractOperator; verbose=false, threshold = default_threshold(A),directsolver=sparseQR_solver, options...)
    if !(directsolver isa typeof(sparseQR_solver))
        error("Only sparse QR possible as direct solver")
    end
    dict = src(A)
    g = grid(dest(A))
    Zt = AZ_Zt(ap; verbose=verbose, threshold=threshold, options...)
    dualdict = dest(Zt)
    AAZA = sparseAAZAmatrix(dict, dualdict, g)
    step1 = directsolver(ArrayOperator(AAZA,src(A),dest(A)); verbose=verbose, threshold=threshold, options...)
end



end
