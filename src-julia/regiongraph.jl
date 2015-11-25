function regiongraph(aff,seg,max_segid)
    """
    # given a segmentation and an affinity graph, create a region graph
    # edge weight is maximum affinity between two regions in the affinity graph
    # edge to self weighted by maximum affinity within a region
    """
    (xdim,ydim,zdim)=size(seg)
    edges=[Dict{UInt32,Float64}() for i=1:max_segid]    # Array of Dict
#    edges=[i => Dict{UInt32,Float64}() for i=1:max_segid]  #Dict of Dict, not sure which is more efficient

    for z=1:zdim
        for y=1:ydim
            for x=1:xdim
                if seg[x,y,z]!=0   # ignore singleton voxels
                    if ( (x > 1) && seg[x-1,y,z]!=0 )
#                        if seg[x,y,z]!=seg[x-1,y,z]
                            (id1,id2) = minmax(seg[x,y,z], seg[x-1,y,z])
                            curr = get(edges[id1],id2,0)
                            edges[id1][id2] = max(curr, aff[x,y,z,1])
#                        end
                    end
                    if ( (y > 1) && seg[x,y-1,z]!=0 )
#                        if seg[x,y,z]!=seg[x,y-1,z] 
                            (id1,id2) = minmax(seg[x,y,z], seg[x,y-1,z])
                            curr = get(edges[id1],id2,0)
                            edges[id1][id2] = max(curr, aff[x,y,z,2])
#                        end
                    end
                    if ( (z > 1) && seg[x,y,z-1]!=0 )
#                        if seg[x,y,z]!=seg[x,y,z-1]
                            (id1,id2) = minmax(seg[x,y,z], seg[x,y,z-1])
                            curr = get(edges[id1],id2,0)
                            edges[id1][id2] = max(curr, aff[x,y,z,3])
#                        end
                    end
                end
            end
        end
    end

    rg = []
    for id1::UInt32 = 1:max_segid
        for (id2, weight) in edges[id1]
            push!(rg, (weight, id1, id2))
        end
    end

    println("Region graph size: ", length(rg))
    sort!(rg,rev=true)  # for stability, should algorithm be specified?
    println("Sorted")

    return rg
end
