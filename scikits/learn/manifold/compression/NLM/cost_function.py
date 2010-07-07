
# Matthieu Brucher
# Last Change : 2007-07-18 14:14

import numpy
import itertools

def dist2hd(x,y):
   """
   Generate a 'coordinate' of the solution at a time
   """
   d = numpy.zeros((x.shape[0],y.shape[0]),dtype=x.dtype)
   for i in xrange(x.shape[1]):
       diff2 = x[:,i,None] - y[:,i]
       diff2 **= 2
       d += diff2
   numpy.sqrt(d,d)
   return d

class CostFunction(object):
  """
  Cost function for the NLM algorithm
  """
  def __init__(self, distances, *args, **kwargs):
    """
    Saves the distances to approximate
    """
    self.distances = distances
    self.factor = numpy.sum(distances)
    self.len = len(self.distances)

  def __call__(self, parameters):
    """
    Computes the cost for a parameter
    """
    params = parameters.reshape((self.len, -1))
    d = dist2hd(params, params)
    d = (d-self.distances)**2/self.distances
    d[numpy.where(numpy.isnan(d))] = 0
    return self.factor * numpy.sum(d)

  def gradient(self, parameters):
    """
    Gradient of this cost function
    """
    params = parameters.reshape((self.len, -1))
    d = dist2hd(params, params)

    grad = numpy.zeros(params.shape)
    for (g, x, d_a, d_r) in itertools.izip(grad, params, d, self.distances):
      temp = 2 * (x - params).T * (d_a-d_r)/(d_r*d_a)
      temp[numpy.where(numpy.isnan(temp))] = 0
      g[:]= numpy.sum(temp, axis=1)
    return grad.ravel()
