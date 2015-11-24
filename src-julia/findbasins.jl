function findbasins(sag::Array{UInt8,3})
    """
    input is the steepest ascent graph
    all paths are unique, except in maximal plateaus
    find a vertex that has not been assigned
    follow path in directed graph
    (1) If you reach a voxel previously assigned a segment ID, 
    then assign voxels on path to same ID.
    (2) If you run out of voxels, assign voxels on path to new ID.
    """
    seg=convert(Array{UInt32,3}, sag)
    counts0 = 0  # number of singleton vertices
    counts=[]   # will store voxel counts for each segment
    bfs=[]
#    counts=Array{Int64,1}(0)     # why is this slower?
#    bfs=Array{Int64,1}(0)

    (xdim,ydim,zdim) = size(seg) 
    const dir = [-1, -xdim, -xdim*ydim, 1, xdim, xdim*ydim]
    const dirmask  = [0x01, 0x02, 0x04, 0x08, 0x10, 0x20]
    
    # definitions should be generalized to UInt64
    # MSB indicates whether voxel has been assigned a segment ID
    const high_bit = 0x80000000::UInt32
    const low_bits = 0x7FFFFFFF::UInt32

    next_id = 1   # initialize segment ID
    for idx in eachindex(seg)
        if seg[idx] == 0   # singleton vertex (no edges at all)
            seg[idx] |= high_bit 
            counts0 += 1;
        elseif (seg[idx] & high_bit)==0  # not assigned
            push!(bfs,idx)     # enqueue
            seg[idx] |= 0x40    # mark as visited
            
            bfs_index = 1  # follow trajectory starting from idx
            while ( bfs_index <= length(bfs) )
                me = bfs[bfs_index]
                for d = 1:6
                    if ( seg[me] & dirmask[d] ) !=0  # outgoing edge
                        him = me + dir[d]  # target of edge
                        if ( seg[him] & high_bit ) !=0 # if already assigned
                            # assign every vertex in queue to same
                            counts[ seg[him] & low_bits ] += length(bfs);
                            for it in bfs
                                seg[it] = seg[him]  # including high bit
                            end
                            bfs = []  # empty queue
                            break
                        elseif ( ( seg[him] & 0x40 ) == 0 )  # not yet visited
                            seg[him] |= 0x40;    # mark as visited
                            push!(bfs,him)    # enqueue
                        end
                    end
                end
                bfs_index += 1      # go to next vertex in queue
            end

            if length(bfs) != 0     # new segment has been created
                push!(counts,length(bfs))
                for it in bfs
                    seg[it] = high_bit | next_id    # assign a segment ID
                end
                next_id += 1
                bfs = []
            end
        end
    end

    println("found: ", (next_id-1)," components\n")

    for idx in eachindex(seg)
        seg[idx] &= low_bits     # clear MSB
    end

    (seg, counts, counts0)
end
