using Watershed

export wsseg, watershed, mergert!, mergert

function watershed(affs, low=0.3, high=0.8, thresholds=[(600,0.3)], dust_size=1000)
    println("watershed, low: $low, high: $high")
    sag = steepestascent(affs, low, high)
    divideplateaus!(sag)
    (seg, counts, counts0) = findbasins(sag)
    rg = regiongraph(affs, seg, length(counts))
    new_rg = mergeregions(seg, rg, counts, thresholds, dust_size)
    rt = mst(new_rg, length(counts))
    return (seg, rt)
end

function mergert(seg, rt, thd=0.5)
    # the returned segmentation
    ret = deepcopy(seg)
    mergert!(ret, rt, thd)
    return ret
end

function mergert!(seg, rt, thd=0.5)
    # get the ralative parent dict
    pd = Dict()
    # initialized as children and parents
    num = 0
    for t in rt
        a, c, p = t
        @assert p>0 && c>0
        if a >= thd
            num = num + 1
        end
        pd[c] = (p, a)
    end
    println("number of merging edges: $num")
    println("total number: $(length(rt))")

    # get the relative root id
    # root dict
    rd = Dict()
    # root set
    rset = IntSet()
    for t in rt
        # get affinity and segment IDs of parent and child
        a, c0, p0 = t
        c = c0
        p = p0
        while a >= thd && haskey(pd,p)
            # next pair of child and parent
            a = pd[p][2]
            p = pd[p][1]
        end
        if p != p0
            rd[c0] = p
            push!(rset, p)
        end
    end
    println("number of trees: $(length(rset))")

    #println("root dict: $rd")
    num = 0
    # set the segment id as relative root id
    for i in eachindex(seg)
        # segment id
        sid = seg[i]
        if haskey(rd, sid)
            #println("root: $(UInt64(root))")
            seg[i] = rd[sid]
            num = num + 1
        end
    end
    println("really merged edges: $num")
end

function wsseg2d(affs, low=0.3, high=0.9, thresholds=[(256,0.3)], dust_size=100, thd_rt=0.5)
    seg = zeros(UInt32, size(affs)[1:3] )
    for z in 1:size(affs,3)
        seg[:,:,z], rt = watershed(affs[:,:,z,:], low, high, thresholds, dust_size)
        seg[:,:,z] = mergert(seg[:,:,z], rt, thd_rt)
    end
    return seg
end

function wsseg(affs, dim = 3, low=0.3, high=0.9, thresholds=[(256,0.3)], dust_size=100, thd_rt=0.5)
    @assert dim==2 || dim==3
    if dim==2
        return wsseg2d(affs, low, high, thresholds, dust_size, thd_rt)
    else
        seg, rt = watershed(affs, low, high, thresholds, dust_size)
        seg = mergert(seg, rt, thd_rt)
        return seg
    end
end
