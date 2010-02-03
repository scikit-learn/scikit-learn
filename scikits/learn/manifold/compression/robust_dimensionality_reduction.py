
"""
Robust optimization with a specific cost function
"""

import numpy
import numpy.linalg
import math

from scikits.optimization import *
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
  def __init__(self, nb_coords):
    self.nb_coords = nb_coords

  def __call__(self, parameters):
    points = parameters.reshape((-1, self.nb_coords))
    means = numpy.mean(points, axis = 0)
    return (points - means).ravel()

class AddNoise(object):
  """
  Adds a small amount of noise to each coordinates
  """
  def __init__(self, nb_coords, function, temp = 1.):
    self.nb_coords = nb_coords
    self.function = function
    self.temp = temp / 10

  def __call__(self, parameters):
    cost = self.function(parameters)
    print cost
    points = parameters.reshape((-1, self.nb_coords))
    cost /= points.shape[0]**2 * points.shape[1] # mean cost per distance
    cost = math.sqrt(cost)
    print cost

    if (cost * self.temp)**(1/8.) > 0:
      points += numpy.random.normal(loc = 0, scale = (cost * self.temp)**(1/8.), size = points.shape)
    self.temp /= 1.5

    return points.ravel()

def optimize_cost_function(distances, function, nb_coords = 2, **kwargs):
  """
  Computes a new coordinates system that respects the distances between each point
  Parameters :
    - distances is the distances to respect
    - nb_coords is the number of remaining coordinates
    - epsilon is a small number
    - sigma is the percentage of distances below which the weight of the cost function is diminished
    - x1 is the percentage of distances '' which the weight of the cost function is diminished
    - x2 is the percentage of distances that indicates the limit when differences between estimated and real distances are too high and when the cost becomes quadratic
  """
  import pickle
  function = function(distances, nb_coords, **kwargs)
  if 'x0' in kwargs:
    x0 = kwargs['x0']
  else:
    x0 = numpy.zeros(distances.shape[0] * nb_coords)#numpy.random.normal(0., math.sqrt(variance), distances.shape[0] * nb_coords)

  err = numpy.seterr(invalid='ignore')

  optimi = optimizer.StandardOptimizerModifying(
    function = function,
    step = step.GradientStep(),
    criterion = criterion.criterion(ftol = 0.00000001, iterations_max = 10000),
    x0 = x0,
    line_search = line_search.FixedLastStepModifier(step_factor = 4., line_search = line_search.FibonacciSectionSearch(alpha_step = 1., min_alpha_step = 0.0001)),
    pre_modifier = AddNoise(nb_coords, function),
    post_modifier = Modifier(nb_coords))

  optimal = optimi.optimize()
  optimal = optimal.reshape(-1, nb_coords)

  numpy.seterr(**err)

  return optimal
