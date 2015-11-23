function mst(rg, max_segid)
    """
    compute maximal spanning tree from weighted graph
    rg - region graph with edges sorted by descending weight
    max_segid - largest node ID in graph
    regiontree - maximal spanning tree
    if edges of rg are sorted by ascending weight, will compute minimal spanning tree 
    """
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
    order = zeros(max_segid)   # will contain numbering of vertices
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
}
