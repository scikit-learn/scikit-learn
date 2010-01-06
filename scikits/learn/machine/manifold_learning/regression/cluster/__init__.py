
"""
Clustering Module
"""

# Matthieu Brucher
# Last Change : 2007-11-08 15:55

from _modified_general_clustering import *

__all__ = ['GeneralCluster']

def test(level=-1, verbosity=1):
  from numpy.testing import NumpyTest
  return NumpyTest().test(level, verbosity)
