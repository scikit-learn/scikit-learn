
"""
Maximum Likelihood Piecewise Linear Mapping Regression module
"""

import math
import numpy
import numpy.linalg as linalg
from numpy.random import shuffle
import random
import copy

import PLMR
import logging

class MLPLMR(PLMR.PLMR):
  """
  Regression with piecewise linear functions
  Uses ML or mean square error (same error for every piecewise function)
  """
  def __init__(self, points, coords, criterion = None, **kwargs):
    """
    Initializes the regression
    - points are the initial points
    - coords are the coordinates that will be used
    - neighbors is the number of neighboor used for determining a plan's equation
    - random_variable is the kid of random variable that will be used for estimation, it is supposed to be identical for every piecewise function
    - criterion is the stopping criterion
    """
    if not criterion:
      from scikits.optimization import criterion
      self.criterion = criterion.ModifiedAICCriterion(-0.00001, 1000, (coords.shape[-1] * points.shape[-1]) / (30 * numpy.std(points)))
    else:
      self.criterion = criterion
    super(MLPLMR, self).__init__(points, coords, **kwargs)
    self.iteration = 0

  def learn(self):
    """
    Tries to learn the model
    """
    self.state = {'old_value' : 1.e3000, 'old_parameters' : [], 'iteration' : 0}
    iterMax = 100

    self.belonging_vector = numpy.zeros(len(self.coords), dtype = numpy.int)
    self.equations = [numpy.array((0,))]
    self.updateEquations()
    self.state['new_parameters'] = copy.deepcopy(self.equations)
    self.state['new_value'] = -self._getLogLikelihood()

    while not self.criterion(self.state):
      candidates = self.getBestCandidates()

      candidate = random.randint(0, len(candidates) - 1)
      oldBV = self.belonging_vector.copy()

      self.findEquationAround(candidates[candidate])
      tempBV = self.belonging_vector.copy()
      self.updateBV()

      underiter = 0
      while (tempBV != self.belonging_vector).any():
        if underiter > iterMax:
          self.belonging_vector = oldBV
          self.updateEquations()
          break
        tempBV = self.belonging_vector.copy()
        self.pruneEquations()
        self.updateEquations()
        self.updateBV()
        underiter += 1

      if (numpy.max(self.belonging_vector) < self.coords.shape[0]/100) and (self.ensure_connexity()):
        while (tempBV != self.belonging_vector).any():
          if underiter > iterMax:
            self.belonging_vector = oldBV
            self.updateEquations()
            break
          tempBV = self.belonging_vector.copy()
          self.pruneEquations()
          self.updateEquations()
          self.updateBV()
          underiter += 1

      self.state['iteration'] +=1
      self.state['old_parameters'] = self.state['new_parameters']
      self.state['old_value'] = self.state['new_value']
      self.state['new_parameters'] = copy.deepcopy(self.equations)
      self.state['new_value'] = -self._getLogLikelihood()
      logging.debug("Equation(s): %d, likelihood: %f", len(self.equations), self.state['new_value'])

    self.belonging_vector = oldBV
    self.updateEquations()

    self.computeError()
    del self.points
    self.RBFF = [self.createRBF(numpy.where(self.belonging_vector == plan)[0]) for plan in range(0, len(self.equations))]

  def updateBV(self):
    """
    Updates the belonging vector
    """
    self.computeError()
    errors = numpy.array([[self.random_variable.get_log_likelihood(point - numpy.dot(coord, equation)) for (coord, point) in zip(self.coords, self.points)] for equation in self.equations])
    self.belonging_vector = numpy.argmax(errors, axis=0)

  def _getLogLikelihood(self):
    """
    Computes the log-likelihood of the model
    """
    errors = self.computeResiduals();
    logErrors = numpy.array(errors)**2
    return -numpy.sum(logErrors)

  def getBestCandidates(self, n = 30):
    """
    Returns the n best group of candidates for a new plan
    """
    errors = self.computeResiduals();
    logErrors = numpy.log(numpy.array(errors)**2)

    candidates = numpy.zeros(len(logErrors))
    for point in range(0, len(errors)):
      candidates[point] = -numpy.sum(logErrors[self.graph[point]])

    return candidates.argsort()[:n]
