function regiongraph{T}(aff::Array{T,4},seg,max_segid)
    # how to constrain so fourth dimension has size three?
    """
    # given a segmentation and an affinity graph, create a region graph
    # edge weight is maximum affinity between two regions
    # self edge weight is maximum affinity within a region
    # output is (weights, vertices)
    """
    (xdim,ydim,zdim)=size(seg)

    # edge list representation
    edges=Dict{Tuple{UInt32,UInt32},T}()
    # keys are vertex pairs (i,j) where i \leq j
    # values are edge weights
    # efficiency is competitive with Array of Dicts and code is simpler
    
    low = convert(T,0)  # choose a value lower than any affinity in the region graph
    
    for z=1:zdim
        for y=1:ydim
            for x=1:xdim
                if seg[x,y,z]!=0   # ignore singleton voxels
                    if ( (x > 1) && seg[x-1,y,z]!=0 )
                        p = minmax(seg[x,y,z], seg[x-1,y,z])
                        edges[p] = max(get(edges,p,low), aff[x,y,z,1])
                    end
                    if ( (y > 1) && seg[x,y-1,z]!=0 )
                        p = minmax(seg[x,y,z], seg[x,y-1,z])
                        edges[p] = max(get(edges,p,low), aff[x,y,z,2])
                    end
                    if ( (z > 1) && seg[x,y,z-1]!=0 )
                        p = minmax(seg[x,y,z], seg[x,y,z-1])
                        edges[p] = max(get(edges,p,low), aff[x,y,z,3])
                    end
                end
            end
        end
    end

    # separate weights and vertices in two arrays
    nedges = length(edges)
    weights = zeros(Float64,nedges)
    vertices = zeros(UInt32,2,nedges)
    i = 1
    for (p, weight) in edges
        weights[i]=weight
        vertices[:,i]=collect(p)
        i +=1
    end
    println("Region graph size: ", nedges)

    # sort both arrays so that weights decrease
    p = sortperm(weights,rev=true)
    weights = weights[p]
    vertices = vertices[:,p]
    return weights, vertices
end
