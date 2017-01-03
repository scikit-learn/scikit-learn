"""Fixture module to skip the unsupervised_learning.rst doctest for 
versions of SciPy earlier than 0.12.0. 
"""
from sklearn.utils.testing import SkipTest
from sklearn.utils.fixes import sp_version

def setup_module(module):
    if sp_version < (0, 12):
        raise SkipTest("Skipping because SciPy version earlier than 0.12.0 and "
                       "thus does not include the scipy.misc.face() image.")
