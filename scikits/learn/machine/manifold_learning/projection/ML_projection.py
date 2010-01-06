
"""
Projection with ML on a piecewise linear function module
"""

from scikits.optimization import *

import numpy
import numpy.linalg as linalg

import scipy.optimize

import logging

__all__ = ['MLProjection']

class ML(object):
  """
  Cost function based on ML probabilities that must be maximized
  """
  def __init__(self, probaEps, Y, equation, mask):
    """
    Constructs the probabibly as a product
    - probaY is the probability on espilon
    - Y is the point to project
    - equation is the matrix to go from X to Y (Y=WX)
    """
    self.probaEps = probaEps
    self.Y = numpy.squeeze(Y)
    self.equation = equation
    self.mask = mask

  def __call__(self, x):
    """
    Computes the proba for the sample x
    """
    return -self.probaEps(self.Y - numpy.dot(x, self.equation), self.mask)

  def gradient(self, x):
    """
    Computes the gradient of a posteriori sample
    """
    grad = numpy.dot(self.probaEps.gradient(self.Y - numpy.dot(x, self.equation), self.mask), self.equation.T)
    return numpy.squeeze(grad)

class MLProjection(object):
  """
  Class that will handle the projection
  - PLMR is an instance of PLMR or that satisfies its attribute interface
  """
  def __init__(self, PLMR):
    self.PLMR = PLMR
    self.PLMRcost = self.PLMR.getLogLikelihood

  def project(self, point, mask=1):
    """
    Project a new point on the manifold described by PLMR
    """
    candidates = {}
    for equation in range(len(self.PLMR.equations)):
      (cost, coord, epsilon, index) = self.projectOnPlan(point, equation, mask)
      candidates[cost] = (coord, epsilon, index)
    c = numpy.array(candidates.keys())
    best = numpy.nanmin(c)
    logging.debug("Likelihood: %f, coordinates: %s, error: %s, projection: %s, original: %s", best, candidates[best][0], epsilon, point-epsilon, point)
    return (candidates[best][0], self.function(candidates[best][0], candidates[best][2]), best)

  def function(self, coord, index):
    """
    Computes the point on the manifold
    """
    return numpy.dot(coord, self.PLMR.equations[index])

  def projectOnPlan(self, point, equation, mask):
    """
    Projects a point on a plan and returns the coordinates on the plan with the reconstruction error
    """
    centered_point = numpy.asmatrix(point - self.PLMR.equations[equation][-1])
    centered_equation = numpy.asmatrix(self.PLMR.equations[equation][:-1])

    coord = numpy.asarray(centered_point * centered_equation.T * linalg.inv(centered_equation * centered_equation.T))
    coord = self.computeBest(coord.squeeze(), numpy.asarray(centered_point), numpy.asarray(centered_equation), self.PLMR.RBFF[equation], mask)
    coordbis = numpy.append(coord.squeeze(), [1])

    cost = self.PLMRcost(coordbis, point, equation = equation)
    reconstruct = numpy.dot(coordbis, self.PLMR.equations[equation])
    epsilon = point - reconstruct
    logging.debug("Likelihood: %f; coordinates: %s", cost, coord)
    return(cost, coordbis, epsilon, equation)

  def computeBest(self, coord, point, equation, RBFF, mask):
    """
    Computes the best coordinates with maximization of the a posteriori probability of X and the error
    """
    function = ML(self.PLMR.random_variable.RBF, point, equation, mask)
    opt = optimizer.StandardOptimizer(function = function, step = step.GradientStep(), criterion = criterion.criterion(ftol = 0.0001, gtol = 0.0001, iterations_max = 200), x0 = numpy.squeeze(coord), line_search = line_search.FibonacciSectionSearch(min_alpha_step=0.0001))
    return opt.optimize()
    #return scipy.optimize.fmin(function, x0 = coord)
