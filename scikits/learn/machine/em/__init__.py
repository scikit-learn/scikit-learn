#! /usr/bin/env python
# Last Change: Sun Jul 22 01:00 PM 2007 J

from info import __doc__

from gauss_mix import GmParamError, GM
from gmm_em import GmmParamError, GMM, EM
from online_em import OnGMM as _OnGMM

__all__ = filter(lambda s:not s.startswith('_'), dir())

from numpy.testing import NumpyTest
test = NumpyTest().test

def test_suite(*args):
    # XXX: this is to avoid recursive call to itself. This is an horrible hack,
    # I have no idea why infinite recursion happens otherwise.
    if len(args) > 0:
        import unittest
        return unittest.TestSuite()
    np = NumpyTest()
    np.testfile_patterns.append(r'test_examples.py')
    return np.test(level = -10, verbosity = 5)
