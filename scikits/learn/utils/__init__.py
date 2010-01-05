#! /usr/bin/env python
# Last Change: Sat Jul 21 02:00 PM 2007 J
from preprocessing import scale, nanscale, Scaler, NanScaler

from numpy.testing import NumpyTest
__all__ = ['scale', 'nanscale', 'Scaler', 'NanScaler']
test = NumpyTest().test

def test_suite(*args):
    # XXX: this is to avoid recursive call to itself. This is an horrible hack,
    # I have no idea why infinite recursion happens otherwise.
    if len(args) > 0:
        import unittest
        return unittest.TestSuite()
    return NumpyTest().test(level = -10)
