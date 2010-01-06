
"""
Projection with MAP on a piecewise linear function module
"""

from scikits.optimization import *
import numpy
import numpy.linalg as linalg
import scipy.optimize

import math

import ML_projection

__all__ = ['MAPProjection']

class APosteriori(object):
  """
  Cost function based on a posteriori probabilities that must be maximized
  """
  def __init__(self, probaX, probaEps, Y, equation, mask):
    """
    Constructs the probaibly as a product
    - probaX is the probability on X
    - probaY is the probability on espilon
    - Y is the point to project
    - equation is the matrix to go from X to Y (Y=WX)
    """
    self.probaX = probaX
    self.probaEps = probaEps
    self.Y = numpy.squeeze(Y)
    self.equation = equation
    self.mask = mask

  def __call__(self, x):
    """
    Computes the proba for the sample x
    """
    return -(self.probaEps(self.Y - numpy.dot(x, self.equation), self.mask) + self.probaX(x))

  def gradient(self, x):
    """
    Computes the gradient of a posteriori sample
    """
    grad = numpy.dot(self.probaEps.gradient(self.Y - numpy.dot(x, self.equation), self.mask), self.equation.T) - self.probaX.gradient(x)
    if numpy.isnan(grad).any():
      grad = numpy.zeros(grad.shape)
    return numpy.squeeze(grad)

class MAPProjection(ML_projection.MLProjection):
  """
  Class that will handle the projection
  - PLMR is an instance of PLMR or that satisfies its attribute interface
  - neighboors is the number of neigboors to use
  """
  def __init__(self, PLMR):
    self.PLMR = PLMR
    self.PLMRcost = self.PLMR.get_MAP

  def computeBest(self, coord, point, equation, RBFF, mask):
    """
    Computes the best coordinates with maximization of the a posteriori probability of X and the error
    """
    function = APosteriori(RBFF, self.PLMR.random_variable.RBF, point, equation, mask)
    opt = optimizer.StandardOptimizer(function = function, step = step.GradientStep(), criterion = criterion.criterion(ftol = 0.0001, gtol = 0.0001, iterations_max = 200), x0 = numpy.squeeze(coord), line_search = line_search.FibonacciSectionSearch(min_alpha_step=.00001))
    return opt.optimize()
    #return scipy.optimize.fmin(function, x0 = coord)
