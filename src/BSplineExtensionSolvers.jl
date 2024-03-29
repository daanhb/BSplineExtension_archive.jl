module BSplineExtensionSolvers

using SparseArrays, StaticArrays, FrameFun.BasisFunctions, FrameFun.BasisFunctions.GridArrays,
    LinearAlgebra, FrameFun.ExtensionFrames, DomainSets
import FrameFun.BasisFunctions: linearized_apply!, apply!, GenericSolverOperator,
    VectorizingSolverOperator
using FrameFun.FrameFunInterface
using FrameFun.BasisFunctions.GridArrays: boundary_mask, ModCartesianIndices

include("nonzero_cols.jl")
include("nonzero_rows.jl")

export BSplineExtensionSolver
"""
    struct BSplineExtensionSolver{T} <: BasisFunctions.VectorizingSolverOperator{T}

This solver is efficient and effective in the first step of the AZ algorithm of a compact dictionary only.
For this dictionary the function `nonzero_coefficients` should be much smaller than the length of the dictionary.

See also:  [`nonzero_coefficients`](@ref)


    BSplineExtensionSolver(M::DictionaryOperator; directsolver=:qr, crop=false, crop_tol=0, verbose=false; options...)

# Arguments
- `M`: The operator to be solved. Its source should be an `ExtensionFrame`
    (with a `support` and a `superdict`) and
    its destination a GridBasis with collocation points in the source support.
    The indices of the nonzero columns of `M` should be equal to the indices of the
    source superdict elements that evaluate non-zero on
    the source support boundary. An example is AZ'A-A where A is `AZ_A` and Z'
    is `AZ_Zt` of the platform `ExtensionFramePlatform(CDBSplinePlatform(), 0.0..0.5)`

## Keywords
- `directsolver::Symbol = :qr`: The direct solver to use to solve the left-over
    truncated Operator `M`. See `FrameFun`.
- `crop::Bool = true`: Truncate the nonzero rows too.
- `crop_tol::Number = 0`: If `crop` is true, truncate the rows with elements smaller than `crop_tol`.
- `sparse::Bool = false`: convert the cropped matrix to a sparse matrix.
- `verbose::Bool = false`: Print method information. See `FrameFun`

# Examples
The easiest use of `BSplineExtensionSolver` is providing it just with (AZ'A-A)
```jldocs
julia> P = ExtensionFramePlatform(CDBSplinePlatform(), 0.0..0.5); N = 30;

julia> plunge = plungeoperator(P,N); A = AZ_A(P,N); Zt = AZ_Zt(P,N);

julia> M = plunge*A;

julia> S = BSplineExtensionSolver(M;crop=true)
BSplineExtensionSolver



julia> b = samplingoperator(P,N)*exp;

julia> x1 = S*plunge*b;

julia> x2 = Zt*(b-A*x1);

julia> F = DictFun(dictionary(P,N), x1 + x2)
DictFun{Float64,Float64}(A 1-dimensional Expansion with 30 degrees of freedom.
Basis: Extension frame
)
```

More often it will be used as a parameter for `FrameFun` functionality. The following
    does the same as the previous example:
```jldocs
julia> P = ExtensionFramePlatform(CDBSplinePlatform(), 0.0..0.5); N = 30;

julia> Fun(exp, P, N; REG = BSplineExtensionSolver)
DictFun{Float64,Float64}(A 1-dimensional Expansion with 30 degrees of freedom.
Basis: Extension frame
)

```
"""
struct BSplineExtensionSolver{T} <: VectorizingSolverOperator{T}
    op      :: DictionaryOperator{T}
    grid_res     :: DictionaryOperator{T}
    dict_ext     :: DictionaryOperator{T}
    sol     :: VectorizingSolverOperator{T}

    dict_scratch
    grid_scratch

    src_linear  ::  Vector{T}
    dest_linear ::  Vector{T}

    function BSplineExtensionSolver(M::DictionaryOperator{T};
            crop=true, crop_tol=0, directsolver=:qr, verbose=false, lazy=false, sparse=false, options...) where T

        _nonzero_cols = nonzero_cols(basis(src(M)), supergrid(grid(dest(M))), support(src(M)))
        dict_resop = IndexRestrictionOperator(src(M), src(M)[_nonzero_cols], _nonzero_cols)

        verbose && println("BSplineExtensionSolver: Restrict columns (coefficients) from $(size(src(M))) to $(length(_nonzero_cols))")

        if sparse
            verbose && println("BSplineExtensionSolver: Sparsifying")
            m = SparseArrays.sparse(M*dict_resop').A
            verbose && println("BSplineExtensionSolver: fill $(nnz(m)/prod(size(m))*100)%")
        else
            s = zeros(src(M))
            d = zeros(dest(M))
            m = Array{eltype(M)}(undef, size(M,1), length(_nonzero_cols))

            for (j,i) in enumerate(_nonzero_cols )
                s[i] = 1
                apply!(M, d, s)
                s[i] = 0
                copyto!(m, size(M,1)*(j-1)+1, d, 1)
            end
        end
        if crop
            verbose && println("BSplineExtensionSolver: Restricting rows...")
            nz_rows = findall(nonzero_rows(m,size(dest(M));nonzero_tol=crop_tol))[:]
            I = LinearIndices(size(dest(M)))[nz_rows]

            if length(nz_rows) < size(M,1)
                m = m[I,:]
                grid_resop = IndexRestrictionOperator(dest(M), GridBasis{coefficienttype(dest(M))}(grid(dest(M))[I]),nz_rows)

                verbose && println("BSplineExtensionSolver: Restrict rows (collocation points) from $(size(dest(M))) to $(length(nz_rows))")
                if sparse
                    verbose && println("BSplineExtensionSolver: fill $(nnz(m)/prod(size(m))*100)%")
                end
            else
                grid_resop = IdentityOperator(dest(M))
                verbose && println("BSplineExtensionSolver: Rows are not restricted")
            end
        else
            grid_resop = IdentityOperator(dest(M))
        end

        new{T}(M, grid_resop, dict_resop',
            lazy ? GenericSolverOperator(ArrayOperator(m), ArrayOperator(m')) :
                FrameFunInterface.directsolver(ArrayOperator(m); directsolver=directsolver, verbose=verbose, options...),
            zeros(length(_nonzero_cols)), zeros(dest(grid_resop)),
            zeros(T, length(dest(M))), zeros(T, length(src(M))))
    end
end

LinearAlgebra.qr(M::SparseMatrixCSC, ::Val{true};opts...) = qr(M;opts...)

export truncated_size
"""
    truncated_size(op::BSplineExtensionSolver)

The size of the smaller system. The one that has to be solved with a direct solver.
"""
truncated_size(op::BSplineExtensionSolver) = size(inv(op.sol))

function linearized_apply!(op::BSplineExtensionSolver, dest::Vector, src::Vector)
    apply!(op.grid_res, op.grid_scratch, src)
    apply!(op.sol, op.dict_scratch, op.grid_scratch)
    apply!(op.dict_ext, dest, op.dict_scratch)
end

export sparseQR_solver
sparseQR_solver(op::DictionaryOperator; options...) =
    GenericSolverOperator(op, sparseqr_factorization(op))

sparseqr_factorization(op::DictionaryOperator) = qr(sparse(op).A)
end
