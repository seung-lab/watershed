# load dependencies
include("steepestascent.jl")
include("divideplateaus.jl")
include("findbasins.jl")
include("regiongraph2.jl")
include("mergeregions2.jl")
include("mst.jl")

using HDF5
aff=h5read("../../out1/out1.affinity.h5","out1")
println("read affinity graph")
low = .3
high = .9
@time sag=steepestascent(aff,low,high);
println("created steepest ascent graph")
@time divideplateaus!(sag);
println("divided plateaus")
@time (seg, counts, counts0)=findbasins(sag);
println("found basins")
@time (w,v)=regiongraph2(aff,seg,length(counts));
println("created region graph")
@time new_rg=mergeregions2(seg,w,v,counts,[(256,.3)]);
println("merged regions")
@time rt=mst(new_rg,length(counts));
println("mst")

