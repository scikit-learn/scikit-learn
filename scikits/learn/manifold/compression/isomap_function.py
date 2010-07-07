
# Matthieu Brucher
# Last Change : 2008-04-07 16:27

import numpy
import itertools

from tools import dist2hd

class CostFunction(object):
  """
  Cost function for the Isomap algorithm
  """
  def __init__(self, distances, *args, **kwargs):
    """
    Saves the distances to approximate
    """
    self.distances = distances
    self.len = len(self.distances)

  def __call__(self, parameters):
    """
    Computes the cost for a parameter
    """
    params = parameters.reshape((self.len, -1))
    d = dist2hd(params, params)
    diff_d = d**2-self.distances**2
    diff_d -= diff_d.mean(axis=0)[:,None]
    diff_d -= diff_d.mean(axis=1)[None,:]
    d = diff_d**2
    return numpy.sum(d)

  def gradient(self, parameters):
    """
    Gradient of this cost function
    """
    params = parameters.reshape((self.len, -1))
    d = dist2hd(params, params)

    grad = numpy.zeros(params.shape)
    for (g, x, d_a, d_r) in itertools.izip(grad, params, d, self.distances):
      temp = 4 * (d_a**2-d_r**2) * (x - params).T
      temp[numpy.where(numpy.isnan(temp))] = 0
      g[:]= numpy.sum(temp, axis=1)
    return grad.ravel()
