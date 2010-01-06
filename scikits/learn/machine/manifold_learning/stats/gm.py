
# Matthieu Brucher
# Last Change : 2008-03-26 14:05

"""
Laplace Estimator
"""

import math
import numpy

import kernels

class IsotropicGemanMcClureVariable(object):
  """
  An isotropic Geman-McClure variable
  """
  def __init__(self):
    """
    Initialization of the RV
    """
    self.samples = None

  def setup(self):
    """
    Create the RBF
    """
    self.RBF = kernels.IsotropicGemanMcClureKernel(self.mean, self.std)

  def get(self):
    """
    Returns a dictionary with the elements needed to create another kernel
    """
    return {'mean':self.mean, 'std':self.std}

  def set(self, d):
    """
    Updates self with a new set of values
    """
    self.__dict__.update(d)

  def clean(self):
    """
    Cleans the samples
    """
    self.samples = None

  def getLogLikelihood(self, x, **kwargs):
    """
    Returns the likelihood of a point
    """
    return self.RBF(x, **kwargs)
