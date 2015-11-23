# load dependencies
include("steepestascent.jl")
include("divideplateaus.jl")
include("findbasins.jl")

function watershed(aff,low,high)
    sag=steepestascent(aff,low,high)
    divideplateaus!(sag)
    (seg, counts, counts0)=findbasins(sag)
    return seg, counts, counts0
end
