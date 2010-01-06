
"""
Robust optimization with a specific cost function
"""

# Matthieu Brucher
# Last Change : 2008-04-07 19:02

import numpy
import numpy.linalg
import math

from toolbox.optimizers import *
import cost_function

class Recorder(object):
  def __init__(self):
    self.elements = []

  def __call__(self, **state):
    self.elements.append(state.copy())
    del self.elements[-1]['function']

class Modifier(object):
  """
  Recenters the points on each axis
  """
  def __init__(self, nbCoords):
    self.nbCoords = nbCoords

  def __call__(self, parameters):
    points = parameters.reshape((-1, self.nbCoords))
    means = numpy.mean(points, axis = 0)
    return (points - means).ravel()

class AddNoise(object):
  """
  Adds a small amount of noise to each coordinates
  """
  def __init__(self, nbCoords, function, temp = 1.):
    self.nbCoords = nbCoords
    self.function = function
    self.temp = temp / 10

  def __call__(self, parameters):
    cost = self.function(parameters)
    print cost
    points = parameters.reshape((-1, self.nbCoords))
    cost /= points.shape[0]**2 * points.shape[1] # mean cost per distance
    cost = math.sqrt(cost)
    print cost

    if (cost * self.temp)**(1/8.) > 0:
      points += numpy.random.normal(loc = 0, scale = (cost * self.temp)**(1/8.), size = points.shape)
    self.temp /= 1.5

    return points.ravel()

def optimize_cost_function(distances, function, nbCoords = 2, **kwargs):
  """
  Computes a new coordinates system that respects the distances between each point
  Parameters :
    - distances is the distances to respect
    - nbCoords is the number of remaining coordinates
    - epsilon is a small number
    - sigma is the percentage of distances below which the weight of the cost function is diminished
    - x1 is the percentage of distances '' which the weight of the cost function is diminished
    - x2 is the percentage of distances that indicates the limit when differences between estimated and real distances are too high and when the cost becomes quadratic
  """
  import pickle
  function = function(distances, nbCoords, **kwargs)
  if 'x0' in kwargs:
    x0 = kwargs['x0']
  else:
    x0 = numpy.zeros(distances.shape[0] * nbCoords)#numpy.random.normal(0., math.sqrt(variance), distances.shape[0] * nbCoords)

  err = numpy.seterr(invalid='ignore')

  optimi = optimizer.StandardOptimizerModifying(
    function = function,
    step = step.GradientStep(),
    criterion = criterion.criterion(ftol = 0.00000001, iterations_max = 10000),
    x0 = x0,
    line_search = line_search.FixedLastStepModifier(step_factor = 4., line_search = line_search.FibonacciSectionSearch(alpha_step = 1., min_alpha_step = 0.0001)),
    pre_modifier = AddNoise(nbCoords, function),
    post_modifier = Modifier(nbCoords))

  optimal = optimi.optimize()
  optimal = optimal.reshape(-1, nbCoords)

  numpy.seterr(**err)

  return optimal
