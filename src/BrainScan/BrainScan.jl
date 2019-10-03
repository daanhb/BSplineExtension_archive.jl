module BrainScan
include("BrainScanGrids.jl")
using .BrainScanGrids


using MAT
function getdata()
    file = matopen(joinpath(splitdir(@__FILE__())[1], "mridata.mat"))
    data = read(file, "mridata")
    close(file)
    data
end

export rhs
function rhs(;normalized=true, opts...)
    data = getdata()
    grid = grid(data;isin=k->k!=0,opts...)
    b = data[subindices(grid)]
    if normalized
        normalize(b)
    else
        b
    end
end

end
