# load dependencies
using StatsBase

export wsseg, watershed, aff2sgm, atomicseg, mergerg!, mergerg, rg2segmentPairs

function _baseseg(aff::AffinityMap,
            low::AbstractFloat,
            high::AbstractFloat,
            thresholds::Vector,
            dust_size::Int; is_threshold_relative=false)
  if is_threshold_relative
    info("use percentage threshold")
    if length(aff) > 3*1024*1024*128
        h = StatsBase.fit(Histogram,
                        aff[1:1024,1:1024,1:128,:][:]; nbins = 1000000)
    else
      h = StatsBase.fit(Histogram, aff[:]; nbins = 1000000)
    end
    low  = _percent2thd(h, low)
    high = _percent2thd(h, high)
    for i = 1:length( thresholds )
      thresholds[i] = tuple(thresholds[i][1], _percent2thd(h, thresholds[i][2]))
    end
  end
  info("absolute watershed threshold: low: $low, high: $high, thresholds: $(thresholds)")
  # this seg is a steepest ascent graph, it was named as such for in-place computation to reduce memory comsuption
  println("steepestascent...")
  seg = steepestascent(aff, low, high)
  println("divideplateaus...")
  divideplateaus!(seg)
  println("findbasins!")
  (seg, counts, counts0) = findbasins!(seg)
  println("regiongraph...")
  rg = regiongraph(aff, seg, length(counts))
  println("mergeregions...")
  new_rg = mergeregions!(seg, rg, counts, thresholds, dust_size)
  return seg, new_rg, counts
end

function atomicseg(aff::AffinityMap,
            low::AbstractFloat=0.1,
            high::AbstractFloat=0.8,
            thresholds::Vector=[(800,0.2)],
            dust_size::Int=600; is_threshold_relative=false)
    seg, rg, counts = _baseseg(aff, low, high, thresholds, dust_size; is_threshold_relative = is_threshold_relative)
    return seg
end

function watershed(aff::AffinityMap, low::AbstractFloat=0.1, high::AbstractFloat=0.8,
                    thresholds::Vector=[(800,0.2)], dust_size::Int=600; is_threshold_relative=false)
    seg, new_rg, counts = _baseseg(aff, low, high, thresholds, dust_size; is_threshold_relative = is_threshold_relative)
    rg = mst(new_rg, length(counts))
    return (seg, rg)
end

function mergerg(seg::Segmentation, rg::RegionGraph, thd::AbstractFloat=0.5)
    # the returned segmentation
    ret = deepcopy(seg)
    mergerg!(ret, rg, thd)
    return ret
end

function mergerg!(seg::Segmentation, rg::RegionGraph, thd::AbstractFloat=0.5)
    # get the ralative parent dict
    pd = Dict()
    # initialized as children and parents
    num = 0
    for t in rg
        a, c, p = t
        @assert p>0 && c>0
        if a >= thd
            num = num + 1
        end
        pd[c] = (p, a)
    end
    println("number of merging edges: $num")
    println("total number: $(length(rg))")

    # get the relative root id
    # root dict
    rd = Dict()
    # root set
    rset = IntSet()
    for t in rg
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

function _wsseg2d(affs, low=0.3, high=0.9, thresholds=[(256,0.3)], dust_size=100, thd_rt=0.5)
    seg = zeros(UInt32, size(affs)[1:3] )
    for z in 1:size(affs,3)
        seg[:,:,z], rg = watershed(affs[:,:,z,:], low, high, thresholds, dust_size)
        seg[:,:,z] = mergerg(seg[:,:,z], rg, thd_rg)
    end
    return seg
end

function wsseg(affs, dim = 3, low=0.3, high=0.9, thresholds=[(256,0.3)], dust_size=100, thd_rg=0.5)
    @assert dim==2 || dim==3
    if dim==2
        return _wsseg2d(affs, low, high, thresholds, dust_size, thd_rg)
    else
        seg, rg = watershed(affs, low, high, thresholds, dust_size)
        seg = mergerg(seg, rg, thd_rg)
        return seg
    end
end

"""
transform rg to dendrogram for omnification
"""
function rg2segmentPairs(rg::RegionGraph)
    N = length(rg)
    segmentPairAffinities = zeros(Float32, N)
    segmentPairs = zeros(UInt32, N,2)

    Threads.@threads for i in 1:N
        segmentPairAffinities[i] = rg[i][1]
        segmentPairs[i,1] = rg[i][2]
        segmentPairs[i,2] = rg[i][3]
    end
    return segmentPairs, segmentPairAffinities
end


#=doc
.. function::
   Transform ralative threshold to absolute threshold
   Args:
   -
=#

function _percent2thd(h::StatsBase.Histogram, rt::AbstractFloat)
  # total number
  tn = sum(h.weights)
  # the rank of voxels corresponding to the threshold
  rank = tn * rt
  # accumulate the voxel number
  avn = 0
  for i in eachindex(h.weights)
      avn += h.weights[i]
      if avn >= rank
          return h.edges[1][i]
      end
  end
end
