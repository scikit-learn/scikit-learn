#!/usr/bin/env python 

help = """Projects new samples on a manifold described by a regression.

Options:
  -h prints this help

Usage:
  projection.py [projectionkind [datafile [regressionfile [projectedfile]]]]

  - projectionkind is one of the projection algorithm in scikits.learn.machine.manifold_learning.projection (MLProjection, MAPProjection, ...)
  - datafile is the data file to project (default = swissroll.samples.pickled)
  - regressionfile is the model file (default = swissroll.regressed.pickled)
  - projectedfile is the output file (default = swissroll.projected.pickled)
"""

import sys
import pickle
import numpy

from scikits.learn.machine.manifold_learning import projection

if len(sys.argv) > 1:
  if sys.argv[1] == "-h":
    print help
    exit()
  projectionkind = sys.argv[1]
else:
  projectionkind = 'MAPProjection'

if len(sys.argv) > 2:
  datafile = sys.argv[2]
else:
  datafile = "swissroll.samples.pickled"

if len(sys.argv) > 3:
  regressionfile = sys.argv[3]
else:
  regressionfile = "swissroll.regressed.pickled"

if len(sys.argv) > 4:
  projectedfile = sys.argv[4]
else:
  projectedfile = "swissroll.projected.pickled"

print "Importing samples dataset %s" % datafile
f = open(datafile)
data = pickle.load(f)

print "Importing model %s" % regressionfile
f = open(regressionfile)
model = pickle.load(f)

print "Projection using %s" % projectionkind
projectionalgo = getattr(projection, projectionkind)
projection_model = projectionalgo(model)

projecteds = numpy.zeros((0, data.shape[1]))
for sample in data:
  (coord, projected, best) = projection_model.project(sample)
  projecteds = numpy.vstack((projecteds, projected[None,:]))

print "Saving results in %s" % projectedfile
f = open(projectedfile, 'w')
pickle.dump(projecteds, f)
