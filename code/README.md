Watershed
=======

What is it?
------------

This code produces an executable to run watershed segmentation with single linkage clustering as outlined in the attached paper.

The latest version
------------

Should be master branch until we start tagging.

Installation
------------

###Required Libraries

|Library|Ubuntu package name|
|-------|-------------------|
|[boost](http://www.boost.org/)|libboost-all-dev|
###Compilation
```
./make.sh
```

Code Documentation
------------

The algorithm is generally split into 4 stages:

1. [Watershed](#watershed)
1. [Region Graph](#region-graph)
1. [Merge Regions](#merge-regions)
1. [Maximal Spanning Tree](#maximal-spanning-tree)

### Watershed
`(./bin/watershed)`
This function takes the Affinity graph and creates the initial watershed segmentation.

####Inputs
#####`xSize (size_t)`
######Dimension in x direction
#####`ySize (size_t)`
######Dimension in y direction
#####`zSize (size_t)`
######Dimension in z direction
#####`aff (float)`
######Affinity graph represented as a 4-Dimensional float array of [xSize][ySize][zSize][3].
#####`lowv (float)`
######Minimum threshold for watershed
#####`highv (float)`
######Maximum threshold for watershed
####Outputs
#####`seg (uint32_t)`
######3-Dimensional uint32_t array [xSize][ySize][zSize] representing the segmentId for each voxel
#####`counts (size_t)`
######Number of voxels for each segmentId, index is the segmentId
      
### Region Graph
`(./regionGraph)`
####Inputs
#####`aff (float)`
######Affinity graph represented as a 4-Dimensional float array of [xSize][ySize][zSize][3].
#####`seg (uint32_t)`
######3-Dimensional uint32_t array [xSize][ySize][zSize] representing the segmentId for each voxel
#####`max_segid (size_t)`
######Maximum segment id value
####Outputs
#####`rg Array<tuple<float, uint32_t, uint32_t>>`
######Watershed basin graph of nodes and weights - sorted and ordered such that maximum weight is at the beginning of the list

### Merge Regions
`(./bin/mergeRegions)`
Apply size dependent single linkage clustering to merge regions in a region graph. Merged regions are replaced a single region in the region graph and segmentation. TODO add note about func

####Inputs
#####`seg (uint32_t)`
######3-Dimensional uint32_t array [xSize][ySize][zSize] representing the segmentId for each voxel
#####`rg Array<tuple<float, uint32_t, uint32_t>>` 
######Region graph as a list of edges: weight, segmentId, segmentId
#####`counts (size_t)`
######Number of voxels for each segmentId, index is the segmentId
#####`func (void*)` *Not yet implemented, uses only a constant above threshold instead*
######A function with an `()` operator that scales the maximum threshold 
#####`thold (size_t)`
######Maximum size threshold to allow merging two regions together
#####`lowt (size_t)`
######Minimum size threshold needed in order to rewrite the old segmentation ids to the new segmentation ids
####Outputs
#####`seg (uint32_t)` 
######3-Dimensional uint32_t array [xSize][ySize][zSize] representing the new segmentIds
#####`rg Array<tuple<float, uint32_t, uint32_t>>`  
######Region graph of nodes and weights

### Maximal Spanning Tree
`(./bin/maximalSpanningTree)`
Given a graph, find the maximal spanning tree.
####Inputs
#####`rg Array<tuple<float, uint32_t, uint32_t>>`
######Graph as a list of edges: weight, nodeId, nodeId.  Sorted by descending weight. 
#####`max_seg_id (size_t)`
######Maximum segment id
####Outputs
#####`mst Array<tuple<float, uint32_t, uint32_t>>`
######Graph as a list of edges: weight, nodeId, nodeId.  Sorted by descending weight. The child node is on the left and the parent node is on the right. Because this is an mst, the child node on the left will be unique.

Usage
-------
All steps can be automatically performed with the `./bin/rws` binary. Inputs are as indicated from `--help` option.
```
$ ./bin/runWatershedFull --help
Basic command to run watershed on affifinity graph with single linkage clustering. Options:

Generic options:
  --help                print help messages

Input options:
  --inputFile arg       Filename of affinity graph input file.
                          Data: Affinity graph values
                          Data Format:  4-D float array of 
                                       [3][xSize][ySize][zSize] *need to 
                                       confirm this* read in 
                                       fortran_storage_order
  --xSize arg           Dimension in X Direction.
  --ySize arg           Dimension in Y Direction.
  --zSize arg           Dimension in Z Direction.

Watershed Options:
  --lowv arg (=0.300000012)             Minimum threshold for watershed.
  --highv arg (=0.899999976)            Maximum threshold for watershed.
  --enableMerge arg (=1)                Enable merge region step for single 
                                        linkage clustering
  --lowt arg (=256)                     Minimum merge size
  --thold arg (=256)                    Maximum merge size (calculated from 
                                        --func
  --funcName arg (=constant)            Merge thresholding function.
                                        ** NOT IMPLEMENTED Defaults to 
                                        const_above_threshold(.3, thold)**
                                        Example inputs:
                                          --funcName=const --funcArg1=.3 
                                        --funcArg2=1000
                                          --funcName=linear --funcArg1=100
                                          --funcName=square --funcArg1=100
                                          --funcName=power --funcArg1=2 
                                        --funcArg2=5000
  --funcArg1 arg (=0.29999999999999999) Argument 1 for merge thresholding 
                                        function
  --funcArg2 arg                        Argument 2 for merge thresholding 
                                        function
  --funcArg3 arg                        Argument 3 for merge thresholding 
                                        function

Watershed Options:
  --outFileSegment arg (=ws.segment.data.out)
                                        Filename of the segmentation output 
                                        file.
                                          Data: Segmentation Ids for each voxel
                                          Data Format: 3-Dimensional uint32_t 
                                                       array [xSize][ySize][zSi
                                                       ze] 
  --outFileDendPairs arg (=ws.dend_pairs)
                                        Filename of the MST Dendrogram pairs 
                                        file. Sorted by weight first as 
                                        specified in the dendValues file. 
                                        Additionally, the child node which is 
                                        unique will be on the left
                                          Data: Graph Edges for MST (weights in
                                                another file)
                                          Data Format: 2-D array of uint32_t 
                                                       pairs i.e. [[SegId1, 
                                                       SegId2][SegId1, 
                                                       SegId3]...]
  --outFileDendValues arg (=ws.dend_values)
                                        Filename of the MST Dendrogram values 
                                        file.
                                          Data: Graph Edge Weight Probabilities
                                          Data Format: 1-D array of decreasing 
                                                       floats corresponding to 
                                                       the pair of SegIds 
                                                       specified in the 
                                                       dendPairs file

Usage Examples:
    runWatershedFull ws.affinity.data 256 256 256 0.3 0.9 250 10 ws.segment.data ws.dend_pairs ws.dend_values
```

