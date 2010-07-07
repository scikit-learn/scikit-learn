
# Matthieu Brucher
# Last Change : 2007-10-29 15:00

"""
A gaussian kernel
"""

import math
import numpy

class IsotropicGaussianKernel(object):
  """
  An isotropic gaussian kernel
  """
  def __init__(self, loc, scale, *args, **kwargs):
    """
    Initializes the kernel:
    - loc is the center of the kernel
    - scale is the size of the kernel (float)
    """
    self.loc = loc
    self.scale = scale

    self.factor = -math.log(2 * numpy.pi * scale**2)*(len(loc) / 2.)

  def __call__(self, x, mask = 1, **kwargs):
    """
    Computes the log-pdf at x for this kernel
    """
    xp = (x-self.loc) / self.scale
    return self.factor - 1 / 2. * numpy.inner(xp, xp)

  def gradient(self, x, mask = 1, **kwargs):
    """
    Computes the gradient of the kernel for x
    """
    return -1/self.scale * (x - self.loc)

class AnisotropicGaussianKernel(object):
  """
  An anisotropic gaussian kernel
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

    self.factor = -1/2. * math.log(2 * numpy.pi)*len(loc) - numpy.sum(numpy.log(scale))

  def __call__(self, x, mask = 1, **kwargs):
    """
    Computes the log-pdf at x for this kernel
    """
    xp = (x-self.loc) / self.scale
    return self.factor - 1/2. * numpy.inner(xp * mask, xp * mask)

  def gradient(self, x, mask = 1, **kwargs):
    """
    Computes the gradient of the kernel for x
    """
    return -1/self.scale * (x - self.loc) * mask
