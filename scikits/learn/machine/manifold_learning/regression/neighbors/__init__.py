
# Matthieu Brucher
# Last Change : 2008-04-15 10:42

"""
Neighbors module
"""

from neighbors import *
from utilities import *

__all__ = ['Neighbors', 'Kneighbors', 'Parzen', 'create_graph']

def test(level=-1, verbosity=1):
  from numpy.testing import NumpyTest
  return NumpyTest().test(level, verbosity)
