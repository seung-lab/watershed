# load dependencies
include("steepestascent.jl")
include("divideplateaus.jl")
include("findbasins.jl")
include("regiongraph.jl")
include("mergeregions.jl")
include("mst.jl")

using HDF5
@time aff=h5read("../../../out1/out1.affinity.h5","out1")
println("read affinity graph")
low = .3
high = .9
@time sag=steepestascent(aff,low,high);
println("created steepest ascent graph")
@time divideplateaus!(sag);
println("divided plateaus")
@time (seg, counts, counts0)=findbasins(sag);
println("found basins")
@time bg=regiongraph(aff,seg,length(counts));
println("created basin graph")
@time rg=mergeregions(seg,bg,counts,[(256,.3)]);
println("merged basins")
@time rt=mst(rg,length(counts));
println("mst")

