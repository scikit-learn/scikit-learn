#!/usr/bin/env python 

help = """Compresses a simple dataset.

Options:
  -h prints this help

Usage:
  compression.py [compressionkind [datafile [compresseddatafile]]]

  - compressionkind is one of the compression in scikits.learn.machine.manifold_learning.compression (isomap, LLE, ...)
  - datafile is the data file to compress (default = swissroll.pickled)
  - compresseddatafile is the output file (default = swissroll.compressed.pickled)
"""

import sys
import pickle
import numpy

from scikits.learn.machine.manifold_learning import compression

if len(sys.argv) > 1:
  if sys.argv[1] == "-h":
    print help
    exit()
  compressionkind = sys.argv[1]
else:
  compressionkind = 'isomap'

if len(sys.argv) > 2:
  datafile = sys.argv[2]
else:
  datafile = "swissroll.pickled"

if len(sys.argv) > 3:
  compresseddatafile = sys.argv[3]
else:
  compresseddatafile = "swissroll.compressed.pickled"

print "Importing dataset %s" % datafile
f = open(datafile)
data = pickle.load(f)

print "Compressing using %s" % compressionkind
compressionalgo = getattr(compression, compressionkind)
coords = compressionalgo(data, nb_coords=2)

print "Saving results in %s" % compresseddatafile
f = open(compresseddatafile, 'w')
pickle.dump(coords, f)
