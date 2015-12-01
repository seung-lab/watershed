using DataStructures

doc"""
`MST` - compute maximal spanning tree from weighted graph

     regiontree = mst(rg,max_segid)

* `rg`: region graph as list of edges, array of (weight,id1,id2)
  tuples.  The edges should be presorted so that weights are in
  descending order.
* `max_segid`: largest ID in region graph
* `regiontree`: *maximal* spanning tree (MST) of region graph as list of edges, array of `(weight,id1,id2)` tuples. The vertices in each edge are ordered so that `id2` is unique across edges. 

The MST effectively represents the dendrogram for single-linkage
clustering.  Each edge `(weight,id1,id2)` in the MST represents an
internal vertex of the dendrogram located at height = `weight`, i.e.,
a merging of two clusters with score `weight`.  The MST contains a bit
more information than the dendrogram, because `id1` and `id2` are the
root IDs of the two clusters, i.e., the maximum weight edge between
the two clusters is between their elements `id1` and `id2`.

The code should work for general graphs.  If edges of `rg` are
presorted by ascending weight, the code will compute the *minimal*
spanning tree rather than the maximal spanning tree.  
"""

function mst(rg, max_segid)
    regiontree = []
    edges=[Set{UInt32}() for i=1:max_segid]    # Array of Sets

    # Kruskal's algorithm
    sets = IntDisjointSets(max_segid)
    for e in rg
        (v1,v2) = e[2:3]
        s1 = find_root(sets,v1)
        s2 = find_root(sets,v2)
        
        if s1 != s2
            push!(regiontree,e)
            union!(sets,s1, s2)

            push!(edges[v1],v2)   # only necessary for ordering
            push!(edges[v2],v1)
        end
    end

    # rest is only necessary for ordering the vertex pairs in each edge
    # bfs (level order) tree traversal
    order = zeros(UInt32,max_segid)   # will contain numbering of vertices
    curr = 1

    for i = 1:max_segid
        if order[i] == 0
            bfs = Deque{Int}()
            push!(bfs,i)
            order[i] = curr
            curr += 1

            while length(bfs)>0
                x = front(bfs)
                shift!(bfs)

                for y in edges[x]
                    if order[y] == 0
                        order[y] = curr
                        curr += 1
                        push!(bfs,y)
                    end
                end
            end
        end
    end

    # order all edges
    for i in 1:length(regiontree)
        e = regiontree[i]
        if order[e[3]] < order[e[2]]
            regiontree[i] = (e[1], e[3], e[2])
        end
    end
        
    return regiontree
end
