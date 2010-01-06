
"""
Quadratic optimization
"""

# Matthieu Brucher
# Last Change : 2007-12-10 17:46

import numpy
import numpy.random
import numpy.linalg
import math

from toolbox.optimizers import *
import cost_function

class Modifier(object):
  """
  Recenters the points on each axis
  """
  def __init__(self, nbCoords, function):
    self.nbCoords = nbCoords
    self.function = function

  def __call__(self, parameters):
    print self.function(parameters)
    points = parameters.reshape((-1, self.nbCoords))
    means = numpy.mean(points, axis = 0)
    return (points - means).ravel()

def optimize_cost_function(distances, function, nbCoords = 2, **kwargs):
  """
  Computes a new coordinates system that respects the distances between each point
  Parameters :
    - distances is the distances to respect
    - nbCoords is the number of remaining coordinates
  """

  function = function(distances, nbCoords, **kwargs)
  std = numpy.std(numpy.sqrt(distances))/20000
  x0 = numpy.random.normal(0., std, distances.shape[0] * nbCoords)

  err = numpy.seterr(invalid='ignore')

  optimi = optimizer.StandardOptimizerModifying(
    function = function,
    step = step.NewtonStep(),
    criterion = criterion.criterion(gtol = 0.000001, ftol = 0.000001, iterations_max = 10000),
    x0 = x0,
    line_search = line_search.SimpleLineSearch(), post_modifier = Modifier(nbCoords, function))

  optimal = optimi.optimize()
  optimal = optimal.reshape(-1, nbCoords)

  numpy.seterr(**err)

  return optimal
