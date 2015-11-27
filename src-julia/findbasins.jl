# MSB indicates whether voxel has been assigned a segment ID
function high_bit(x::Type{UInt32})
    return 0x80000000::UInt32
end

function high_bit(x::Type{UInt64})
    return 0x8000000000000000LL::UInt64
end
    
function low_bits(x::Type{UInt32})
    return 0x7FFFFFFF::UInt32
end

function low_bits(x::Type{UInt64})
    return 0x7FFFFFFFFFFFFFFFLL::UInt64
end

function findbasins{T}(sag::Array{T,3})
    seg = copy(sag)
    (seg, counts, counts0) = findbasins!(seg)
    return (seg, counts, counts0)
end

# in-place version
function findbasins!{T}(seg::Array{T,3})   # should restrict T to UInt32, UInt64
    """
    input is the steepest ascent graph
    all paths are unique, except in maximal plateaus
    find a vertex that has not been assigned
    follow path in directed graph
    (1) If you reach a voxel previously assigned a segment ID, 
    then assign voxels on path to same ID.
    (2) If you run out of voxels, assign voxels on path to new ID.
    """
    counts0 = 0  # number of singleton vertices

    (xdim,ydim,zdim) = size(seg) 
    const dir = [-1, -xdim, -xdim*ydim, 1, xdim, xdim*ydim]
    const dirmask  = [0x01, 0x02, 0x04, 0x08, 0x10, 0x20]

    counts = UInt32[]  # will store voxel counts for each segment
    bfs = UInt32[]

    next_id = 1   # initialize segment ID
    for idx in eachindex(seg)
        if seg[idx] == 0   # singleton vertex (no edges at all)
            seg[idx] |= high_bit(T) 
            counts0 += 1;
        elseif (seg[idx] & high_bit(T))==0  # not assigned
            push!(bfs,idx)     # enqueue
            seg[idx] |= 0x40    # mark as visited
            
            bfs_index = 1  # follow trajectory starting from idx
            while ( bfs_index <= length(bfs) )
                me = bfs[bfs_index]
                for d = 1:6
                    if ( seg[me] & dirmask[d] ) !=0  # outgoing edge
                        him = me + dir[d]  # target of edge
                        if ( seg[him] & high_bit(T) ) !=0 # if already assigned
                            # assign every vertex in queue to same
                            counts[ seg[him] & low_bits(T) ] += length(bfs);
                            for it in bfs
                                seg[it] = seg[him]  # including high bit
                            end
                            bfs = UInt32[]  # empty queue
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
                    seg[it] = high_bit(T) | next_id    # assign a segment ID
                end
                next_id += 1
                bfs = UInt32[]
            end
        end
    end

    println("found: ", (next_id-1)," components")

    for idx in eachindex(seg)
        seg[idx] &= low_bits(T)     # clear MSB
    end

    (seg, counts, counts0)
end
