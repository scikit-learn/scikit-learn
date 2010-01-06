
# Matthieu Brucher
# Last Change : 2008-04-07 18:57

"""
Allows to compute the nearest neighbors
"""

import numpy
from tools import dist2hd

def parzen(samples, window_size, **kwargs):
  """
  Creates a list of the nearest neighbors in a Parzen window
  """
  l = []

  d = dist2hd(samples, samples)

  for dist in d:
    wi = numpy.where(dist < neighbors)[0]
    l.append(wi)

  return l

def kneigh(samples, neighbors, **kwargs):
  """
  Creates a list of the nearest neighbors in a K-neighborhood
  """
  l = []

  d = dist2hd(samples, samples)

  for dist in d:
    indices = numpy.argsort(dist)
    l.append(indices[:neighbors])

  return l

def NumpyFloyd(dists):
  """
  Implementation with Numpy vector operations
  """
  for indice1 in xrange(len(dists)):
    for indice2 in xrange(len(dists)):
      dists[indice2, :] = numpy.minimum(dists[indice2, :], dists[indice2, indice1] + dists[indice1, :])
