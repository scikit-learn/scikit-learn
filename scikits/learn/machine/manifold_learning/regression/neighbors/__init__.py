
# Matthieu Brucher
# Last Change : 2008-04-07 14:06

"""
Neighboors module
"""

from neighboors import *
from utilities import *

__all__ = ['Neighboors', 'KNeighboors', 'Parzen', 'create_graph']

def test(level=-1, verbosity=1):
  from numpy.testing import NumpyTest
  return NumpyTest().test(level, verbosity)
