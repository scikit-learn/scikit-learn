#!/usr/bin/env python 

help = """Estimates a regression between a simple dataset.

Options:
  -h prints this help

Usage:
  regression.py [regressionkind [datafile [compresseddatafile [regressionfile]]]]

  - regressionkind is one of the regression estimation in scikits.learn.machine.manifold_learning.regression (PLMR, MLPLMR)
  - datafile is the data file to regress (default = swissroll.pickled)
  - compresseddatafile is the compressed file (default = swissroll.isomap.pickled)
  - regressionfile is the output file (default = swissroll.regressed.pickled)
"""

import sys
import pickle
import numpy

from scikits.learn.machine.manifold_learning import regression

if len(sys.argv) > 1:
  if sys.argv[1] == "-h":
    print help
    exit()
  regressionkind = sys.argv[1]
else:
  regressionkind = 'MLPLMR'

if len(sys.argv) > 2:
  datafile = sys.argv[2]
else:
  datafile = "swissroll.pickled"

if len(sys.argv) > 3:
  compresseddatafile = sys.argv[3]
else:
  compresseddatafile = "swissroll.isomap.pickled"

if len(sys.argv) > 4:
  regressionfile = sys.argv[4]
else:
  regressionfile = "swissroll.regressed.pickled"

print "Importing dataset %s" % datafile
f = open(datafile)
data = pickle.load(f)

print "Importing compressed dataset %s" % compresseddatafile
f = open(compresseddatafile)
coords = pickle.load(f)

print "Regression using %s" % regressionkind
regressionalgo = getattr(regression, regressionkind)
model = regressionalgo(data, coords, neighbors = 9)
model.learn()

print "Saving results in %s" % regressionfile
f = open(regressionfile, 'w')
pickle.dump(model, f)
