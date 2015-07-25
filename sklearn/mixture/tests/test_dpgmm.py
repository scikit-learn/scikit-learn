import unittest
import sys

import numpy as np

from sklearn.mixture import DPGMM, VBGMM
from sklearn.mixture.dpgmm import log_normalize
from sklearn.datasets import make_blobs
from sklearn.utils.testing import assert_array_less, assert_equal
from sklearn.mixture.tests.test_gmm import GMMTester
from sklearn.externals.six.moves import cStringIO as StringIO

np.seterr(all='warn')


def test_class_weights():
    # check that the class weights are updated
    # simple 3 cluster dataset
    X, y = make_blobs(random_state=1)
    for Model in [DPGMM, VBGMM]:
        dpgmm = Model(n_components=10, random_state=1, alpha=20, n_iter=50)
        dpgmm.fit(X)
        # get indices of components that are used:
        indices = np.unique(dpgmm.predict(X))
        active = np.zeros(10, dtype=np.bool)
        active[indices] = True
        # used components are important
        assert_array_less(.1, dpgmm.weights_[active])
        # others are not
        assert_array_less(dpgmm.weights_[~active], .05)


def test_verbose_boolean():
    # checks that the output for the verbose output is the same
    # for the flag values '1' and 'True'
    # simple 3 cluster dataset
    X, y = make_blobs(random_state=1)
    for Model in [DPGMM, VBGMM]:
        dpgmm_bool = Model(n_components=10, random_state=1, alpha=20,
                           n_iter=50, verbose=True)
        dpgmm_int = Model(n_components=10, random_state=1, alpha=20,
                          n_iter=50, verbose=1)

        old_stdout = sys.stdout
        sys.stdout = StringIO()
        try:
            # generate output with the boolean flag
            dpgmm_bool.fit(X)
            verbose_output = sys.stdout
            verbose_output.seek(0)
            bool_output = verbose_output.readline()
            # generate output with the int flag
            dpgmm_int.fit(X)
            verbose_output = sys.stdout
            verbose_output.seek(0)
            int_output = verbose_output.readline()
            assert_equal(bool_output, int_output)
        finally:
            sys.stdout = old_stdout


def test_verbose_first_level():
    # simple 3 cluster dataset
    X, y = make_blobs(random_state=1)
    for Model in [DPGMM, VBGMM]:
        dpgmm = Model(n_components=10, random_state=1, alpha=20, n_iter=50,
                      verbose=1)

        old_stdout = sys.stdout
        sys.stdout = StringIO()
        try:
            dpgmm.fit(X)
        finally:
            sys.stdout = old_stdout


def test_verbose_second_level():
    # simple 3 cluster dataset
    X, y = make_blobs(random_state=1)
    for Model in [DPGMM, VBGMM]:
        dpgmm = Model(n_components=10, random_state=1, alpha=20, n_iter=50,
                      verbose=2)

        old_stdout = sys.stdout
        sys.stdout = StringIO()
        try:
            dpgmm.fit(X)
        finally:
            sys.stdout = old_stdout


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
        _, z = g.score_samples(train_obs)
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
        _, z = g.score_samples(train_obs)
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
