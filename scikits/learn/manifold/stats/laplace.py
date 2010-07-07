
# Matthieu Brucher
# Last Change : 2007-10-30 16:50

"""
Laplace Estimator
"""

import math
import numpy

import kernels

class IsotropicLaplaceVariable(object):
  """
  An isotropic Laplace variable
  """
  def __init__(self):
    """
    Initialization of the RV
    """
    self.samples = None

  def addSample(self, x):
    """
    Add one or more samples to the variable
    """
    if self.samples:
      self.samples = numpy.append(self.samples, x, axis=0)
    else:
      self.samples = x

  def compute(self):
    """
    Computes the parameters of the variable
    """
    self.mean = numpy.mean(self.samples, axis=0)
    self.std = numpy.mean(numpy.abs(self.samples - self.mean)+1e-30)

    self.setup()

  def setup(self):
    """
    Create the RBF
    """
    self.RBF = kernels.IsotropicLaplaceKernel(self.mean, self.std)

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


class AnisotropicLaplaceVariable(object):
  """
  An anisotropic Laplace variable
  """
  def __init__(self):
    """
    Initialization of the RV
    """
    self.samples = None

  def addSample(self, x):
    """
    Add one or more samples to the variable
    """
    if self.samples:
      self.samples = numpy.append(self.samples, x, axis=0)
    else:
      self.samples = x

  def compute(self):
    """
    Computes the parameters of the variable
    """
    self.mean = numpy.mean(self.samples, axis=0)
    self.std = numpy.mean(numpy.abs(self.samples - self.mean)+1e-10, axis=0)

    self.setup()

  def setup(self):
    """
    Create the RBF
    """
    self.RBF = kernels.AnisotropicLaplaceKernel(self.mean, self.std)

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
