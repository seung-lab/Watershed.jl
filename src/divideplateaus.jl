# export divideplateaus!

"""
`DIVIDEPLATEAUS!` - Divide plateaus in steepest ascent graph

     divideplateaus!(sag)

* `sag`: steepest ascent graph (directed and unweighted). `sag[x,y,z]` contains 6-bit number encoding edges outgoing from (x,y,z)

Modify steepest ascent graph so as to

1. Divide non-maximal plateaus into paths that exit as quickly as possible
2. Break ties between multiple outgoing edges

Note this is an in-place modification of `sag`
"""

function divideplateaus!(sag::Array{UInt32,3})
    (xdim,ydim,zdim) = size(sag)
    const dir = Vector{Int64}([-1, -xdim, -xdim*ydim, 1, xdim, xdim*ydim])
    const dirmask  = [0x00001, 0x00002, 0x00004, 0x00008, 0x00010, 0x00020]
    const idirmask = [0x00008, 0x00010, 0x00020, 0x00001, 0x00002, 0x00004]

    # queue all vertices for which a purely outgoing edge exists
    bfs = Int64[]
    sizehint!(bfs, length(sag))
    for idx in 1:Int64(length(sag))
        for d=1:6
            if (sag[idx] & dirmask[d]) != 0x00000000   # outgoing edge exists
                if (sag[idx+dir[d]] & idirmask[d]) == 0x00000000  # no incoming edge
                    sag[idx] |= 0x00040
                    append!(bfs,[idx])
                    break
                end
            end
        end
    end
    gc()
    # divide plateaus
    bfs_index = Int64(1);
    while bfs_index <= length(bfs)
        idx = bfs[bfs_index]
        to_set = 0x00000
        for d=1:6
            if (sag[idx] & dirmask[d]) != 0x00000000    # outgoing edge exists
                if (sag[idx+dir[d]] & idirmask[d]) != 0x00000000  # incoming edge
                    if ( sag[idx+dir[d]] & 0x00040 ) == 0x00000000
                        push!(bfs,idx+dir[d])
                        sag[idx+dir[d]] |= 0x00040
                    end
                else  # purely outgoing edge
                    to_set = dirmask[d];
                end
            end
        end
        sag[idx] = to_set    # picks unique outgoing edge, unsets 0x40 bit
        bfs_index += 1
    end
    # manually release the memory of bfs
    bfs = nothing
    gc()
    return sag
end
