import unittest

import nose

import numpy as np

from scikits.learn.mixture import DPGMM, VBGMM
from scikits.learn.mixture.dpgmm import log_normalize
from .test_gmm import GMMTester



def test_log_normalize():
    v = np.array([0.1, 0.8, 0.01, 0.09])
    a = np.log(2*v)
    assert np.allclose(v, log_normalize(a))


class DPGMMTester(GMMTester):
    model = DPGMM

    def score(self, g, train_obs):
        return g.lower_bound()


class TestDPGMMWithSphericalCovars(unittest.TestCase, DPGMMTester):
    cvtype = 'spherical'


class TestDPGMMWithDiagCovars(unittest.TestCase, DPGMMTester):
    cvtype = 'diag'


class TestDPGMMWithTiedCovars(unittest.TestCase, DPGMMTester):
    cvtype = 'tied'


class TestDPGMMWithFullCovars(unittest.TestCase, DPGMMTester):
    cvtype = 'full'


class VBGMMTester(GMMTester):
    model = VBGMM

    def score(self, g, train_obs):
        return g.lower_bound()


class TestVBGMMWithSphericalCovars(unittest.TestCase, VBGMMTester):
    cvtype = 'spherical'


class TestVBGMMWithDiagCovars(unittest.TestCase, VBGMMTester):
    cvtype = 'diag'


class TestVBGMMWithTiedCovars(unittest.TestCase, VBGMMTester):
    cvtype = 'tied'


class TestVBGMMWithFullCovars(unittest.TestCase, VBGMMTester):
    cvtype = 'full'


if __name__ == '__main__':
    nose.runmodule()
