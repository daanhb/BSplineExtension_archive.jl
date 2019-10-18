using GridArrays.ModCartesianIndicesBase: ModCartesianIndices, nbindexlist

function braininside(data)
    res = falses(size(data))
    for i in 1:size(data, 3)
        @show i
        res[:,:,i] .= inside(repeatedfilterbycount(data[:,:,i] .> .07, .3*ones(8), [10,10,10,3,3,3,1,1]))
    end
    res
end


function inside(lessnoisymask::BitArray{N}) where N
    mask = falses(size(lessnoisymask))
    mask[1] = true

    for i in CartesianIndices(size(mask))
        if !lessnoisymask[i]
            for j in nbindexlist(i, size(lessnoisymask))
                if mask[j]
                    mask[i] = true
                end
            end
        end
    end
    mask[end] = true

    for i in reverse(CartesianIndices(size(mask));dims=1)
        if !lessnoisymask[i]
            for j in nbindexlist(i, size(lessnoisymask))
                if mask[j]
                    mask[i] = true
                end
            end
        end
    end
    .!mask
end

function repeatedfilterbycount(noisymask::BitArray{N}, r::Vector{<:Real}, m::Vector{Int}) where N
    for (ri, mi) in zip(r, m)
        noisymask = filterbycount(noisymask, ri, mi)
    end
    noisymask
end

function filterbycount(noisymask::BitArray{N}, r, m::Int=5) where N
    c = Matrix{Int}(undef, size(noisymask)...)
    for ci in CartesianIndices(size(c))
        if noisymask[ci]
            cnt = 0
            mci = ModCartesianIndices(size(c),CartesianIndex(ci.I.-m),CartesianIndex(ci.I.+m))
            for i in mci
                if noisymask[i]
                    cnt += 1
                end
            end
            c[ci] = cnt
        else
            c[ci] = 0
        end
    end
    c .> r*(2m+1)^N
end


D = BrainScan.BrainScanGrids.braindata();D=D/norm(D,Inf)
DD = zeros(260,260,56)
DD[1:256,1:256,1:54] .= D
m = braininside(DD)
mm = Array(m[1:256,1:256,1:54])
# save("gridmask.jld","mask",Array(m[1:256,1:256,1:54]))

using BSON
bson("test",Dict(:mask=>m[1:256,1:256,1:54]))
