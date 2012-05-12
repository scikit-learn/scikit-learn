import unittest

import nose

import numpy as np

from sklearn.mixture import DPGMM, VBGMM
from sklearn.mixture.dpgmm import log_normalize
from .test_gmm import GMMTester

np.seterr(all='warn')


def test_log_normalize():
    v = np.array([0.1, 0.8, 0.01, 0.09])
    a = np.log(2 * v)
    assert np.allclose(v, log_normalize(a), rtol=0.01)


def do_model(self, **kwds):
    return VBGMM(verbose=False, **kwds)


class DPGMMTester(GMMTester):
    model = DPGMM
    do_test_eval = False

    def score(self, g, train_obs):
        _, z = g.eval(train_obs)
        return g.lower_bound(train_obs, z)


class TestDPGMMWithSphericalCovars(unittest.TestCase, DPGMMTester):
    covariance_type = 'spherical'
    setUp = GMMTester._setUp


class TestDPGMMWithDiagCovars(unittest.TestCase, DPGMMTester):
    covariance_type = 'diag'
    setUp = GMMTester._setUp


class TestDPGMMWithTiedCovars(unittest.TestCase, DPGMMTester):
    covariance_type = 'tied'
    setUp = GMMTester._setUp


class TestDPGMMWithFullCovars(unittest.TestCase, DPGMMTester):
    covariance_type = 'full'
    setUp = GMMTester._setUp


class VBGMMTester(GMMTester):
    model = do_model
    do_test_eval = False

    def score(self, g, train_obs):
        _, z = g.eval(train_obs)
        return g.lower_bound(train_obs, z)


class TestVBGMMWithSphericalCovars(unittest.TestCase, VBGMMTester):
    covariance_type = 'spherical'
    setUp = GMMTester._setUp


class TestVBGMMWithDiagCovars(unittest.TestCase, VBGMMTester):
    covariance_type = 'diag'
    setUp = GMMTester._setUp


class TestVBGMMWithTiedCovars(unittest.TestCase, VBGMMTester):
    covariance_type = 'tied'
    setUp = GMMTester._setUp


class TestVBGMMWithFullCovars(unittest.TestCase, VBGMMTester):
    covariance_type = 'full'
    setUp = GMMTester._setUp


if __name__ == '__main__':
    nose.runmodule()
