module SparseBandedMatrices

using FrameFun.BasisFunctions
using FrameFun.BasisFunctions: VerticalBandedMatrix, HorizontalBandedMatrix, ArrayOperator,
        IndexMatrix, ExtensionIndexMatrix, RestrictionIndexMatrix, OperatorSum
using SparseArrays

import SparseArrays: SparseMatrixCSC, sparse

sparse(A::ArrayOperator) = ArrayOperator(sparse(A.A), src(A), dest(A))
function sparse(A::CompositeOperator)
    ArrayOperator(*(map(x->sparse(x).A,reverse(elements(A)))...), src(A), dest(A))
end
sparse(A::OperatorSum) = ArrayOperator(A.val1*sparse(A.op1).A+A.val2*sparse(A.op2).A, src(A), dest(A))
sparse(A::TensorProductOperator) = ArrayOperator(kron(map(x->sparse(x).A, elements(A))...), src(A), dest(A))

sparse(A::ScalingOperator) = ArrayOperator(sparse(A.A, size(A)...), src(A), dest(A))


colptr(::Type{Ti}, M::VerticalBandedMatrix) where Ti = Ti[1+k*length(M.array) for k in 0:size(M,2)]

function rowval(::Type{Ti},M::VerticalBandedMatrix) where Ti
    r = Vector{Ti}(undef,length(M.array)*size(M,2))
    s = Vector{Ti}(undef,length(M.array))
    offset = Ti(M.offset)
    shift = Ti(M.step)
    @inbounds for j in Ti(0):Ti(size(M,2)-1)
        for i in Ti(0):Ti(length(s)-1)
            s[i+1] = mod(offset+i,size(M,1))+1
        end
        sort!(s)
        copyto!(r, 1+j*length(s), s, 1, length(s))
        offset += shift
    end
    r
end

function nzval(::Type{Ti},::Type{Tv}, M::VerticalBandedMatrix) where {Ti,Tv}
    r = Vector{Tv}(undef,length(M.array)*size(M,2))
    s = Vector{Ti}(undef,length(M.array))
    offset = Ti(M.offset)
    shift = Ti(M.step)
    for j in Ti(0):Ti(size(M,2)-1)
        for i in Ti(0):Ti(length(s)-1)
            s[i+1] = mod(offset+i,size(M,1))+1
        end
        sort!(s)
        for (i,si) in enumerate(s)
            r[j*length(s)+i] = M[si,j+1]
        end
        offset += shift
    end
    r
end


function SparseMatrixCSC{Tv,Ti}(M::VerticalBandedMatrix) where {Tv,Ti}
    @assert !SparseArrays.has_offset_axes(M)
    SparseMatrixCSC(size(M)..., colptr(Ti,M), rowval(Ti,M), nzval(Ti,Tv,M))
end

SparseMatrixCSC{Tv,Ti}(M::HorizontalBandedMatrix) where {Tv,Ti} =
    SparseMatrixCSC(M')'


colptr(::Type{Ti},M::ExtensionIndexMatrix) where {Ti}= collect(1:size(M,2)+1)
rowval(::Type{Ti},M::ExtensionIndexMatrix{T,1}) where {T,Ti}= collect(subindices(M))
rowval(::Type{Ti},M::ExtensionIndexMatrix{T,N}) where {T,Ti,N} = LinearIndices(M.original_size)[subindices(M)]
nzval(::Type{Ti},::Type{Vi},M::IndexMatrix) where {Ti,Vi} = ones(Vi,length(subindices(M)))

function SparseMatrixCSC{Tv,Ti}(M::ExtensionIndexMatrix) where {Tv,Ti}
    @assert !SparseArrays.has_offset_axes(M)
    SparseMatrixCSC(size(M)..., colptr(Ti,M), rowval(Ti,M), nzval(Ti,Tv,M))
end

SparseMatrixCSC{Tv,Ti}(M::RestrictionIndexMatrix) where {Tv,Ti} =
    SparseMatrixCSC(M')'

end
