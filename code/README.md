Watershed
=======

# What is it?

This code produces an executable to run watershed segmentation with single linkage clustering as outlined in the attached paper.

# The latest version

Should be master branch until we start tagging.

# Installation
## Required Libraries
|Library|Ubuntu package name|
|[boost](http://www.boost.org/)|libboost-all-dev|

## Compilation
```
./make.sh
```

# Usage
The algorithm is generally split into 4 stages:
1. [Watershed Segmentation](#Watershed Segmentation)
1. Generate Watershed Basin Graphs
1. Merge Basin for Region Graph Generation
1. Generate Minimum Spanning Tree

## Watershed Sementation
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

