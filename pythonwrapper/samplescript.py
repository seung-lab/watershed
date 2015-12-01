import numpy as np
import PyWatershed
import h5py
import time
import os
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
intermediateOut = 'intermediate.h5'
if os.path.isfile(intermediateOut):
    os.remove(intermediateOut)
f2 = h5py.File(intermediateOut)
f2.create_dataset('main', data = segVolume)
f2.close()
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
print 'done'
print 'elapsed time to run all of watershed: ', time.time() - beginTime
