function regiongraph{T}(aff::Array{T,4},seg,max_segid)
    """
    # given a segmentation and an affinity graph, create a region graph
    # edge weight is maximum affinity between two regions
    # self edge weight is maximum affinity within a region
    # edgelist - output
    """
    (xdim,ydim,zdim)=size(seg)

    # adjacency list representation
    adjacency=Dict{UInt32,T}[Dict{UInt32,T}() for i=1:max_segid]    # Array of Dict
    low = convert(T,0)  # choose a value lower than any affinity in the region graph
    
    for z=1:zdim
        for y=1:ydim
            for x=1:xdim
                if seg[x,y,z]!=0   # ignore singleton voxels
                    if ( (x > 1) && seg[x-1,y,z]!=0 )
                            (id1,id2) = minmax(seg[x,y,z], seg[x-1,y,z])
                            curr = get(adjacency[id1],id2,low)
                            adjacency[id1][id2] = max(curr, aff[x,y,z,1])
                    end
                    if ( (y > 1) && seg[x,y-1,z]!=0 )
                            (id1,id2) = minmax(seg[x,y,z], seg[x,y-1,z])
                            curr = get(adjacency[id1],id2,low)
                            adjacency[id1][id2] = max(curr, aff[x,y,z,2])
                    end
                    if ( (z > 1) && seg[x,y,z-1]!=0 )
                            (id1,id2) = minmax(seg[x,y,z], seg[x,y,z-1])
                            curr = get(adjacency[id1],id2,low)
                            adjacency[id1][id2] = max(curr, aff[x,y,z,3])
                    end
                end
            end
        end
    end

    # convert to edge list representation
    edgelist = Tuple{T,UInt32,UInt32}[]
    for id1::UInt32 = 1:max_segid
        for (id2, weight) in adjacency[id1]
           push!(edgelist, (weight, id1, id2))
        end
    end
    println("Region graph size: ", length(edgelist))

    edgelist = sort!(edgelist,rev=true,lt=isless2)   # sort is slow
    println("Sorted")

    return edgelist
end

# this comparison function gives substantial speedup
function isless2(x::Tuple{Float64,UInt32,UInt32},y::Tuple{Float64,UInt32,UInt32})
    return isless(x[1],y[1])
end
