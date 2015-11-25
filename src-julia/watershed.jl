# load dependencies
include("steepestascent.jl")
include("divideplateaus.jl")
include("findbasins.jl")
include("regiongraph.jl")
include("mergeregions.jl")
include("mst.jl")

function watershed(aff,low,high)
    sag=steepestascent(aff,low,high)
    divideplateaus!(sag)
    (seg, counts, counts0)=findbasins(sag)

    max_segid = length(counts)
    rg = regiongraph(aff,seg, max_segid)

    thresholds=[(200, 0.0)] 
    dust_size = 0
    rg = merge_segments(seg, rg, counts, thresholds, dust_size)
    tree = mst(rg, max_segid)

    # println(tree)
    return seg, tree, max_segid
end
