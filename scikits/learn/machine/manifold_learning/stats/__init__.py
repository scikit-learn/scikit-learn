
# Matthieu Brucher
# Last Change : 2008-03-26 14:06

"""
Stats module
"""

import kernels

from gaussian import *
from laplace import *
from gm import *
from radial_basis_functions_field import *

__all__ = ['IsotropicGaussianVariable', 'AnisotropicGaussianVariable', 'IsotropicLaplaceVariable', 'AnisotropicLaplaceVariable', 'IsotropicGemanMcClureVariable', 'RBFField']

def test(level=1, verbosity=1):
  from numpy.testing import NumpyTest
  return NumpyTest().test(level, verbosity)
