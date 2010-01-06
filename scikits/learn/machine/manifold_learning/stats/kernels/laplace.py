
# Matthieu Brucher
# Last Change : 2007-09-11 14:18

"""
A Laplace kernel
"""

import math
import numpy

class IsotropicLaplaceKernel(object):
  """
  An isotropic Laplace kernel
  """
  def __init__(self, loc, scale, *args, **kwargs):
    """
    Initializes the kernel:
    - loc is the center of the kernel
    - scale is the size of the kernel (float)
    """
    self.loc = loc
    self.scale = scale
    self.invscale = (1/self.scale)

    self.factor = -math.log(2 * scale)*len(loc)

  def __call__(self, x, mask = 1, **kwargs):
    """
    Computes the log-pdf at x for this kernel
    """
    xp = (x-self.loc) * self.invscale
    return self.factor - numpy.sum(numpy.abs(xp))

  def gradient(self, x, mask = 1, **kwargs):
    """
    Computes the gradient of the kernel for x
    """
    return - numpy.sign(x - self.loc) * self.invscale

class AnisotropicLaplaceKernel(object):
  """
  An anisotropic Laplace kernel
  """
  def __init__(self, loc, scale, *args, **kwargs):
    """
    Initializes the kernel:
    - loc is the center of the kernel
    - scale is the size of the kernel (same size as loc)
    """
    assert(len(loc)==len(scale))
    self.loc = loc
    self.scale = scale
    self.invscale = (1/self.scale)

    self.factor = - numpy.log(2 * scale)

  def __call__(self, x, mask = 1, **kwargs):
    """
    Computes the log-pdf at x for this kernel
    """
    xp = (x-self.loc) * self.invscale
    return numpy.sum(self.factor * mask) - numpy.sum(numpy.abs(xp * mask))

  def gradient(self, x, mask = 1, **kwargs):
    """
    Computes the gradient of the kernel for x
    """
    return - numpy.sign(x - self.loc) * self.invscale * mask
