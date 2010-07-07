
# Matthieu Brucher
# Last Change : 2008-03-26 14:53

"""
A Laplace kernel
"""

import math
import numpy

class IsotropicGemanMcClureKernel(object):
  """
  An isotropic Geman-McClure kernel
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

  def __call__(self, x, mask = 1, **kwargs):
    """
    Computes the log-pdf at x for this kernel
    """
    xp = ((x-self.loc) * self.invscale)**2
    return - numpy.sum(xp / (self.scale + xp))

  def gradient(self, x, mask = 1, **kwargs):
    """
    Computes the gradient of the kernel for x
    """
    xp = ((x-self.loc) * self.invscale)
    return - 2 * xp / (self.scale + xp**2)**2
