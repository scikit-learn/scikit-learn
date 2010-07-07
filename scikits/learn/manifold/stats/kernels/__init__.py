
# Matthieu Brucher
# Last Change : 2008-03-26 14:15

"""
Kernels module
"""

from gaussian import *
from laplace import *
from gm import *

__all__ = ['IsotropicGaussianKernel', 'AnisotropicGaussianKernel', 'IsotropicLaplaceKernel', 'AnisotropicLaplaceKernel', 'IsotropicGemanMcClureKernel']

def test(level=-1, verbosity=1):
  from numpy.testing import NumpyTest
  return NumpyTest().test(level, verbosity)
