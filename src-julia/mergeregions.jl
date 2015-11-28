using DataStructures

function mergeregions(seg, weights, vertices, counts, thresholds, dust_size = 0)
    """
    seg - segmentation.  IDs of foreground segments are 1:length(counts).  ID of background is 0 (modified in place)
    rg - region graph.  IDs should be same as in seg, except no zeros
    new_rg - new region graph (output)
    counts - sizes of regions in `seg` (modified in place)
    thresholds - sequence of (size_th,weight_th) pairs to be used for merging
    dust_size - after merging, tiny regions less than dust_size to be eliminated by changing them to background voxels
    """
    sets = IntDisjointSets(length(counts))
    for (size_th,weight_th) in thresholds
        for (weight,id1,id2) in zip(weights,vertices[1,:],vertices[2,:])
            s1 = find_root(sets,id1)
            s2 = find_root(sets,id2)
            if (weight > weight_th) && (s1 != s2)
                if ( (counts[s1] < size_th) || (counts[s2] < size_th) )
                    counts[s1] += counts[s2]
                    counts[s2]  = 0
                    union!(sets,s1,s2)
                    s = find_root(sets,s1)   # this is either s1 or s2
                    (counts[s], counts[s1]) = (counts[s1], counts[s])
                end
            end
        end
    end
    println("Done with merging")

    # define mapping from parents to new segment IDs
    # and apply to redefine counts
    remaps = zeros(UInt32,length(counts))     # generalize to include UInt64
    next_id = 1
    for id = 1:length(counts)
        s = find_root(sets,id)
        if ( (remaps[s] == 0) && (counts[s] >= dust_size) )
            remaps[s] = next_id
            counts[next_id] = counts[s]    # exercise: prove that next_id < counts
            next_id += 1
        end
    end
    resize!(counts,next_id-1)

    # apply remapping to voxels in seg
    # note that dust regions will get assigned to background
    for idx in eachindex(seg)
        if seg[idx] !=0    # only foreground voxels
            seg[idx] = remaps[find_root(sets,seg[idx])]
        end
    end
    println("Done with remapping, total: ", (next_id-1), " regions")

    # apply remapping to region graph
    in_rg = [Set{UInt32}() for i=1:next_id-1]
    new_rg = Array{Tuple{Float64,UInt32,UInt32},1}(0)
    for (weight, id1, id2) in zip(weights,vertices[1,:],vertices[2,:])
        s1 = remaps[find_root(sets,id1)]
        s2 = remaps[find_root(sets,id2)]
        if ( s1 != s2 && s1 !=0 && s2 !=0)  # ignore dust regions
            (s1,s2) = minmax(s1,s2)
            if ~in(s2,in_rg[s1])
                push!(new_rg,(weight, s1, s2))
                push!(in_rg[s1],s2)
            end
        end
    end

    println("Done with updating the region graph, size: ", length(new_rg))
    return new_rg
end

