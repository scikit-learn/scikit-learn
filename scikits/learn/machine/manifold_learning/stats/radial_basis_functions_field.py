
"""
Class for using RBF fields
An RBF is only a kernel function
"""

import math
import numpy

from scikits.optimization import helpers
from scipy.stats import gaussian_kde

class RBFField(helpers.ForwardFiniteDifferences, gaussian_kde):
  """
  A simple RBF field
  """
  def __init__(self, samples, weight, *args, **kwargs):
    """
    Populates the field with RBF for each sample :
    - samples is a sequence of points in the field
    - RBF_type is the type of RBF to use
    - variance is the variance of a RBF in the field
    - args are the formals arguments for the RBF
    - kwargs is the set of arguments of the RBF constructor
    """
    gaussian_kde.__init__(self, dataset = samples.T)
    helpers.ForwardFiniteDifferences.__init__(self, *args, **kwargs)
    self.weight = weight

  def __call__(self, x):
    """
    Computes the log-probability of a new point x
    """
    p = self.weight * self.evaluate(x)
    if p == 0:
      return -1e10000
    return math.log(p)
