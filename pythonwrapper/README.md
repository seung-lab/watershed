Watershed Python wrapper
=======

What is it?
------------

Build a python extension to call cpp watershed from python

The latest version
------------

Should be *dev* branch until we start tagging.

Installation
------------

###Required Libraries
|Library|Ubuntu package name|Notes|
|-------|-------------------|-----|
|[boost](http://www.boost.org/)|libboost-all-dev||
|[Boost.NumPy](https://github.com/ndarray/Boost.NumPy/)|**N/A**|see installation notes below|

To install Boost.NumPy, try doing:

1. ```sudo apt-get install scons```
1. Clone Boost.NumPy into your directory ```git clone https://github.com/ndarray/Boost.NumPy/```
1. Initialize submodule ```SConsChecks``` folder by running ```git submodule update --init --recursive``` from the Boost.NumPy folder
1. Install using ```scons install```
  * Note 1: if you have a funny boost installation directory you will have to add the flag ```--with-boost```
  * Note 2: default installs Boost.NumPy into ```/usr/local``` which is not automatically picked up by g++ without manually pointing to it with compiler flags. if you add ```--prefix=/usr/``` the compiler should automaticallly pick it up without linking extra directories. (though you will still need ```-lboost_numpy```


###Compilation
```
./make.sh
```
* if it fails, might be because the include directories are not correct, check your installed library locations

Sample Usage
------------

```python
import numpy as np
import PyWatershed
import h5py
import time
import os
# dummy data
#dim = [3, 4, 2, 3]
# data = np.reshape(np.arange(1, np.prod(dim)+1, dtype="float32"), dim)

# Example loading from h5 source
filename = '/usr/people/ww12/seungmount/research/kisuklee/Sharing/Jonathan/Piriform/7nmDeeper/VeryDeep2HR/binary_double/affinitybinaries/out1/ws.affinity.h5';
f = h5py.File(filename, 'r')
print 'reading in data'
data = f['out1'][:].astype('float32')

beginTime = time.time()

# Example for running each step individually (slower for python conversions)
print 'done reading in data', data.shape, 'running watershed'
[segVolume, segCounts] = PyWatershed.watershed(data, .3, .9)
print 'running region graph'
regionGraph = PyWatershed.regionGraph(data, segVolume, len(segCounts)-1)
print 'running merge segments'
print 'before merge we have # segcounts', len(segCounts)
print 'before merge we have # of segments', len(np.unique(segVolume))
print 'before merge we have # of edges', len(regionGraph)
PyWatershed.mergeSegments(segVolume, regionGraph, segCounts, 256, .3, 256);
print 'the new segcounts', len(segCounts)
print 'the new segmentation volume has # of labels', len(np.unique(segVolume))
print 'the new region graph has # of edges', len(regionGraph)
print 'running MST'
mergeTree = PyWatershed.mergeTree(regionGraph, len(segCounts)-1)

# Example for running the all steps in watershed with default arguments
config = {'lowv':.3}
[segVolume, mergeTree] = PyWatershed.watershedFull(data, config)

print 'done'
print 'elapsed time to run all of watershed: ', time.time() - beginTime
````
