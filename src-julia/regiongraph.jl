doc"""
`REGIONGRAPH` - create region graph by finding maximum affinity between each pair of regions in segmentation

     rg = regiongraph(aff,seg,max_segid)

* `rg`: region graph as list of edges, array of (weight,id1,id2) tuples. The edges are sorted so that weights are in descending order.
* `aff`: affinity graph (undirected and weighted). 4D array of affinities, where last dimension is of size 3
* `seg`: segmentation.  Each element of the 3D array contains a *segment ID*, a nonnegative integer ranging from 0 to `max_segid`
* `max_segid`: number of segments

The vertices of the region graph are regions in the segmentation.  An
edge of the region graph corresponds to a pair of regions in the
segmentation that are connected by an edge in the affinity graph.  The
weight of an edge in the region graph is the maximum weight of the
edges in the affinity graph connecting the two regions.

The region graph includes every edge between a region and itself.  
The weight of a self-edge is the maximum affinity within the region.

Background voxels (those with ID=0) are ignored.
"""

function regiongraph{Taff,Tseg}(aff::Array{Taff,4},seg::Array{Tseg,3},max_segid)
    (xdim,ydim,zdim) = size(seg)
    @assert size(aff) == (xdim,ydim,zdim,3)

    # edge list representation
    edges = Dict{Tuple{Tseg,Tseg},Taff}()
    # keys are vertex pairs (i,j) where i \leq j
    # values are edge weights
    # efficiency is competitive with Array of Dicts and code is simpler
    
    low = convert(Taff,0)  # choose a value lower than any affinity in the region graph
    
    for z=1:zdim
        for y=1:ydim
            for x=1:xdim
                if seg[x,y,z]!=0   # ignore background voxels
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

    # sorting array of tuples is now efficient in julia 0.5
    rg = Tuple{Taff,Tseg,Tseg}[]
    for (p, weight) in edges
        push!(rg, (weight, p[1], p[2]))
    end
    sort!(rg, by=edge->edge[1], rev=true)  # "by" arg needed for efficiency
    return rg

    # get rid of clunky workaround code
    # separate weights and vertices in two arrays
#    nedges = length(edges)
#    weights = zeros(Taff,nedges)
#    vertices = zeros(Tseg,2,nedges)
#    i = 1
#    for (p, weight) in edges
#        weights[i]=weight
#        vertices[:,i]=collect(p)
#        i +=1
#    end
#    println("Region graph size: ", nedges)

    # sort both arrays so that weights decrease
#    p = sortperm(weights,rev=true)
#    weights = weights[p]
#    vertices = vertices[:,p]

    # repackage in array of typles
#    rg = Array{Tuple{Taff,Tseg,Tseg},1}(nedges)
#    for i = 1:nedges
#        rg[i]= (weights[i], vertices[1,i], vertices[2,i])
#    end

#    return rg
end
