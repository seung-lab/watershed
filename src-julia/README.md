Hierarchical watershed segmentation
=======

1. Segment an affinity graph into watershed basins.  
   This defines a steepest ascent dynamics on the affinity graph, and then finds the basins of attraction. The basins are in one-to-one-correspondence with regional maxima of the affinity graph.
2. Merge basins via size-dependent single linkage clustering to create regions.  
Watershed typically results in severe oversegmentation. The watershed basins are merged via single linkage clustering to create larger regions.
3. Return region graph and its maximal spanning tree (hierarchical segmentation).  
Regions are vertices of the region graph. The weight of the edge between two regions is defined as the maximum weight of the edges between the two regions in the affinity graph.

`aff` - affinity graph associated with 3D grayscale image.  
    size(aff) = (xdim,ydim,zdim, 3)  
    aff[x,y,z,1] is affinity of voxels at [x-1,y,z] and [x,y,z]  
    aff[x,y,z,2] is affinity of voxels at [x,y-1,z] and [x,y,z]  
    aff[x,y,z,3] is affinity of voxels at [x,y,z-1] and [x,y,z]  
`low` - low threshold  
`high` - high threshold  
`seg` - segmentation represented as 3D indexed image  
    size(seg)=(xdim,ydim,zdim)  
`rg` - region graph  
`rt` - maximal spanning tree
       
Julia translation of Zlateski's C++ code for watershed with size-dependent single linkage clustering.


Code Documentation
------------

The algorithm is generally split into 6 stages:

1. [Steepest ascent graph]
1. [Divide plateaus]
1. [Find basins]
1. [Region Graph](#region-graph)
1. [Merge Regions](#merge-regions)
1. [Maximal Spanning Tree](#maximal-spanning-tree)
