# The number of nonzero values in A end A-AZ^*A using
# N = (128, 128, 27)
# and the degrees as showed the table below

#             | nnz(A)    | nnz(AAZA)
# p = (1,1,1) |   9,000,916 |   2,419,821
# p = (2,2,2) |  28,351,071 |  35,717,131
# p = (3,3,3) |  57,603,712 | 150,645,741

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
        normalize(b, Inf)
    else
        b
    end
end

include("BrainPlatforms.jl")
@reexport using .BrainPlatforms

end
