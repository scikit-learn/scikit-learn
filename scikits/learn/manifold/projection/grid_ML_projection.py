
"""
Projection with ML on a piecewise linear function module with a grid
"""

# Matthieu Brucher
# Last Change : 2008-06-11 09:24

import numpy

import scipy.optimize

__all__ = ['GridMLProjection']

class GridMLProjection(object):
  """
  Class that will handle the projection
  - PLMR is an instance of PLMR or that satisfies its attribute interface
  """
  def __init__(self, PLMR):
    self.PLMR = PLMR
    self.mins = numpy.min(PLMR.coords[:,:-1], axis=0)
    self.maxs = numpy.max(PLMR.coords[:,:-1], axis=0)

    self.extremas = tuple([slice(min - (max - min) / 2.,max + (max - min) / 2.,(max - min) / 10.)  for min, max in zip(self.mins, self.maxs)])
    self.PLMRcost = self.PLMR.get_log_likelihood

  def project(self, point, mask=1):
    """
    Project a new point on the manifold described by PLMR
    """
    candidates = {}

    grid = numpy.array(numpy.mgrid[self.extremas]).reshape(len(self.extremas), -1).T

    for coord in grid:
      coords = numpy.append(coord, [1])
      cost = self.PLMRcost(coords, point, mask)
      reconstruct = self.PLMR.get_point(coords)
      epsilon = point - reconstruct
      candidates[cost] = (coords, epsilon, -1)
    c = numpy.array(candidates.keys())
    indices = numpy.argsort(c[numpy.isreal(c)])

    for indice in indices[-5*len(self.mins):]:
      coords = self.optimize(candidates[c[indice]][0][:-1], point, mask)
      cost = self.PLMRcost(coords, point, mask)
      reconstruct = self.PLMR.get_point(coords)
      epsilon = point - reconstruct
      candidates[cost] = (coords, epsilon, -1)
    c = numpy.array(candidates.keys())
    print len(c)
    best = numpy.nanmin(c)
    print best, candidates[best][0]

    return (candidates[best][0], self.PLMR.get_point(candidates[best][0]), best)

  def optimize(self, coords, point, mask):
    """
    Optimizes a set of coordinates
    """
    coord = scipy.optimize.fmin(self.opt_func, coords, [point, mask], disp=0)
    coords = numpy.append(coord, [1])
    return coords

  def opt_func(self, coords, point, mask):
    coord = numpy.append(coords, [1])
    return self.PLMRcost(coord, point, mask)
