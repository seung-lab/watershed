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

Required Libraries
------------

|Library|Ubuntu package name|
|-------|-------------------|
|[boost](http://www.boost.org/)|libboost-all-dev|
## Compilation
```
./make.sh
```

Code Documentation
------------

The algorithm is generally split into 4 stages:

1. [Watershed Segmentation](#watershed-segmentation)
1. [Generate Watershed Basin Graphs](#generate-watershed-basin-graphs)
1. [Generate Region Graph](#generate-region-graph)
1. [Generate Minimum Spanning Tree](#generate-minimum-spanning-tree)

### Watershed Segmentation
This function takes the Affinity graph and creates the initial watershed segmentation.

####Inputs
#####`xSize (size_t)`
######Dimension in x direction
#####`ySize (size_t)`
######Dimension in y direction
#####`zSize (size_t)`
######Dimension in z direction
#####`aff (float)`
######4-Dimensional float array of [xSize][ySize][zSize][3] *TODO need to confirm this* read in fortran_storage_order
#####`lowv (float)`
######Minimum threshold for watershed
#####`highv (float)`
######Maxmimum threshold for watershed
####Outputs
#####`seg (uint32_t)`
######3-Dimensional uint32_t array [xSize][ySize][zSize] representing the segmentId for each voxel
#####`counts (size_t)`
######Number of voxels for each segmentId *TODO confirm/ask later- size of the vector? is the index the segmentId?*
      
### Generate Watershed Basin Graphs
####Inputs
#####`aff (float)`
######4-Dimensional float array of [xSize][ySize][zSize][3] *TODO need to confirm this* read in fortran_storage_order
#####`seg (uint32_t)`
######3-Dimensional uint32_t array [xSize][ySize][zSize] representing the segmentId for each voxel
#####`max_segid (size_t)`
######Maximum segment id value
####Outputs
#####`rg Array<tuple<float, uint32_t, uint32_t>>`
######Watershed basin graph of nodes and weights - sorted and ordered such that maximum weight is at the beginning of the list

### Generate Region Graph
Use the watershed basin graphs and perform single linkage clustering.  Rewrites the original segmentation and region graphs *TODO does it actually do this?*
####Inputs
#####`seg (uint32_t)`
######3-Dimensional uint32_t array [xSize][ySize][zSize] representing the segmentId for each voxel
#####`counts (size_t)`
######Number of voxels for each segmentId *confirm/ask later- size of the vector? is the index the segmentId?*
#####`func (void*)`
######A function with an `()` operator that scales the maximum threshold
#####`thold (size_t)`
######Maximum size threshold to allow merging two regions together
#####`lowt (size_t)`
######Minimum size threshold needed in order to rewrite the old segmentation ids to the new segmentation ids
####Outputs
#####`seg (uint32_t)` *INPUT IS MODIFIED*
######3-Dimensional uint32_t array [xSize][ySize][zSize] representing the new segmentIds
#####`rg Array<tuple<float, uint32_t, uint32_t>>`  *INPUT IS MODIFIED*
######Region graph of nodes and weights *TODO is sort order preserved?*

### Generate Minimum Spanning Tree
####Inputs
#####`rg Array<tuple<float, uint32_t, uint32_t>>`
######Region graph of nodes and weights
#####`max_seg_id (size_t)`
######Maximum segment id

Usage
-------
All steps can be automatically performed with the `./bin/rws` binary. Inputs are as indicated from `--help` option.
```
./bin/rws --help
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
  --xSize arg           Length of X dimension of input file.
  --ySize arg           Length of Y dimension of input file.
  --zSize arg           Length of Z dimension of input file.

Watershed Options:
  --lowv arg (=0.300000012)  Minimum watershed.
  --highv arg (=0.899999976) Maximum thresholding for waterhsed.
  --lowt arg (=256)          Minimum merge size
  --thold arg (=256)         Maximum merge size (calculated from --func
  --funcName arg (=constant) Merge thresholding function.
                             Example inputs:
                               --funcName=const --funcArg1=.3 --funcArg2=1000
                               --funcName=linear --funcArg1=100
                               --funcName=square --funcArg1=100
                               --funcName=power --funcArg1=2 --funcArg2=5000
  --funcArg1 arg             Argument 1 for merge thresholding function
  --funcArg2 arg             Argument 2 for merge thresholding function
  --funcArg3 arg             Argument 3 for merge thresholding function

Watershed Options:
  --outFileSegment arg (=ws.segment.data.out)
                                        Filename of the segmentation output 
                                        file.
                                          Data: Segmentation Ids
                                          Data Format: 3-D float array of 
                                                       [zSize][ySize][xSize] in
                                                       row-major order
  --outFileDendPairs arg (=ws.dend_pairs)
                                        Filename of the MST Dendrogram pairs 
                                        file. Sorted by weight first as 
                                        specified in the dendValues file.
                                          Data: Segmentation Ids
                                          Data Format: 2-D array of uint32_t 
                                                       pairs i.e. [[SegId1, 
                                                       SegId2][SegId1, 
                                                       SegId3]...]
  --outFileDendValues arg (=ws.dend_values)
                                        Filename of the MST Dendrogram values 
                                        file.
                                          Data: Edge Weight Probabilities
                                          Data Format: 1-D array of decreasing 
                                                       floats corresponding to 
                                                       the pair of SegIds 
                                                       specified in the 
                                                       dendPairs file

Usage Examples:
    ws ws.affinity.data 256 256 256 0.3 0.9 250 10 ws.segment.data ws.dend_pairs ws.dend_values
```
See arguments specified from `--help`

