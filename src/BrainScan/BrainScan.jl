# The number of nonzero values in A end A-AZ^*A using
# N = (128, 128, 27)
# and the degrees as showed the table below

#             | nnz(A)    | nnz(AAZA)
# p = (1,1,1) |  26271216 |  12888778
# p = (2,2,2) |  82752678 | 191614327
# p = (3,3,3) | 168133824 | 840293625

# Only for the lowest degree, it is more efficient to use the sparse AZ algorithm

module BrainScan
include("BrainScanGrids.jl")
using .BrainScanGrids, GridArrays, LinearAlgebra, Reexport
using .BrainScanGrids: braindata

export brainrhs
function brainrhs(;normalized=true, opts...)
    data = braindata()
    b = data[subindices(braingrid(;opts...))]
    if normalized
        normalize(b)
    else
        b
    end
end

include("BrainPlatforms.jl")
@reexport using .BrainPlatforms

end
