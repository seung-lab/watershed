Hierarchical watershed segmentation
=======
Segmentation of an "affinity graph" into "watershed basins," followed by single linkage clustering.  The resulting hierarchy can be represented by a "dendrogram."  In this tree, leaves represent watershed basins, and internal vertices represent mergings of clusters. The height of each vertex is the affinity at which the merging occurs. The tendency to produce oversegmentation (basins that are excessively numerous and/or small) is moderated in two ways:
1. Collapsing all subtrees above a height (*high threshold*) to single leaves. The resulting dendrogram defines a hierarchy on segments that are unions of watershed basins.
2. Collapsing subtrees to eliminate leaves below a *size threshold* (but only for vertices above some height).

This implementation is for an affinity graph associated with a 3D image. Each vertex of the graph is an image voxel, and each edge between nearest neighbor pairs of voxels. A voxel has 6 nearest neighbors in the x,y,z,-x,-y, and -z directions, and the graph is said to have 6-connectivity.  The weight of an edge represents the “affinity” of two voxels for each other.  High affinity voxels tend to end up in the same segment, and low affinity voxels in different segments. 

Watershed basins are associated with local *maxima* of the graph, following the sign convention from graph partitioning in which edges between segments are minimized (and therefore edges within segments are maximized). This is potentially confusing, as the sign convention is the opposite from the original watershed definition in image processing.  In that field, watershed basins are associated with local minima, in accord with the metaphor of a "drop of water flowing downhill."

The affinity graph is represented by three images, because the number of edges in the affinity graph is three times the number of voxels in the image (neglecting boundary effects).  We use convolutional networks to generate the affinity graph from the image [Turaga et al. 2010], but other algorithms may be used.  The simplest case is an affinity graph in which the weight of each edge is the maximum of its two voxel values. Then the watershed basins are associated with local maxima of the image, and our watershed on graphs reduces to the conventional watershed on images (except for the flipped sign).

In image processing, many watershed implementations identify the ridgelines that separate basins, and assign the label "0" to voxels on ridgelines.  Our watershed does not do this, because the ridgeline between two basins is regarded as being located on edges rather than voxels. However, our watershed can be made to label some voxels as background as follows.  Edges of the graph with affinity below a *low threshold* are removed. After this operation, some voxels become singletons, completely disconnected from the rest of the graph. These background voxels are given a label of "0" to distinguish them from foreground regions.

The algorithm
-------

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


Functions
------------

The algorithm is generally split into 6 stages:

1. [Steepest ascent graph]
1. [Divide plateaus]
1. [Find basins]
1. [Region Graph](#region-graph)
1. [Merge Regions](#merge-regions)
1. [Maximal Spanning Tree](#maximal-spanning-tree)

Usage
------------
