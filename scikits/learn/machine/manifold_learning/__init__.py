
"""
Manifold Learning Module
"""

# Matthieu Brucher
# Last Change : 2007-06-13 17:40

import compression
import projection
import regression

def test(level=-1, verbosity=1):
  from numpy.testing import NumpyTest
  return NumpyTest().test(level, verbosity)
