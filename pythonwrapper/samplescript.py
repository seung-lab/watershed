import numpy as np
import PyWatershed
import h5py
import time
dim = [3, 4, 2, 3]
# data = np.reshape(np.arange(1, np.prod(dim)+1, dtype="float32"), dim)
# data = np.reshape(np.arange(1, np.prod(dim)+1, dtype="float32"), dim)

filename = '/usr/people/ww12/seungmount/research/kisuklee/Sharing/Jonathan/Piriform/7nmDeeper/VeryDeep2HR/binary_double/affinitybinaries/out1/ws.affinity.h5';
f = h5py.File(filename, 'r')
print 'reading in data'
data = f['out1'][:].astype('float32')

beginTime = time.time()
print 'done reading in data', data.shape, 'running watershed'
[segVolume, segCounts] = PyWatershed.watershed(data, .3, .9)
print 'running region graph'
regionGraph = PyWatershed.regionGraph(data, segVolume, segCounts.shape[0]-1)
print 'running merge segments'
print 'before merge we have # of elements', len(regionGraph)
PyWatershed.mergeSegments(segVolume, regionGraph, segCounts, 256, .3, 256);
print 'the new region graph has # of elements', len(regionGraph)
print 'running MST'
mergeTree = PyWatershed.mergeTree(regionGraph, segCounts.shape[0]-1)
print 'done'
print 'elapsed time to run all of watershed: ', time.time() - beginTime
