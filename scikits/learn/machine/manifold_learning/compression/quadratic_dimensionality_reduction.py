
"""
Quadratic optimization
"""

import numpy
import numpy.random
import numpy.linalg

from scikits.optimizer import StandardOptimizerModifying
import cost_function

class Modifier(object):
  """
  Recenters the points on each axis
  """
  def __init__(self, nb_coords, function):
    self.nb_coords = nb_coords
    self.function = function

  def __call__(self, parameters):
    print self.function(parameters)
    points = parameters.reshape((-1, self.nb_coords))
    means = numpy.mean(points, axis = 0)
    return (points - means).ravel()

def optimize_cost_function(distances, function, nb_coords = 2, **kwargs):
  """
  Computes a new coordinates system that respects the distances between each point
  Parameters :
    - distances is the distances to respect
    - nb_coords is the number of remaining coordinates
  """

  function = function(distances, nb_coords, **kwargs)
  std = numpy.std(numpy.sqrt(distances))/20000
  x0 = numpy.random.normal(0., std, distances.shape[0] * nb_coords)

  err = numpy.seterr(invalid='ignore')

  optimi = StandardOptimizerModifying(
    function = function,
    step = step.NewtonStep(),
    criterion = criterion.criterion(gtol = 0.000001, ftol = 0.000001, iterations_max = 10000),
    x0 = x0,
    line_search = line_search.SimpleLineSearch(), post_modifier = Modifier(nb_coords, function))

  optimal = optimi.optimize()
  optimal = optimal.reshape(-1, nb_coords)

  numpy.seterr(**err)

  return optimal
