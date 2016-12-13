# Important note for the deprecation cleaning of 0.20 :
# All the function and classes of this file have been deprecated in 0.18.
# When you remove this file please also remove the related files
# - 'sklearn/mixture/dpgmm.py'
# - 'sklearn/mixture/gmm.py'
# - 'sklearn/mixture/test_gmm.py'
import unittest
import sys

import numpy as np

from sklearn.mixture import DPGMM, VBGMM
from sklearn.mixture.dpgmm import log_normalize
from sklearn.datasets import make_blobs
from sklearn.utils.testing import assert_array_less, assert_equal
from sklearn.utils.testing import assert_warns_message, ignore_warnings
from sklearn.mixture.tests.test_gmm import GMMTester
from sklearn.externals.six.moves import cStringIO as StringIO
from sklearn.mixture.dpgmm import digamma, gammaln
from sklearn.mixture.dpgmm import wishart_log_det, wishart_logz


np.seterr(all='warn')


@ignore_warnings(category=DeprecationWarning)
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


@ignore_warnings(category=DeprecationWarning)
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


@ignore_warnings(category=DeprecationWarning)
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


@ignore_warnings(category=DeprecationWarning)
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


@ignore_warnings(category=DeprecationWarning)
def test_digamma():
    assert_warns_message(DeprecationWarning, "The function digamma is"
                         " deprecated in 0.18 and will be removed in 0.20. "
                         "Use scipy.special.digamma instead.", digamma, 3)


@ignore_warnings(category=DeprecationWarning)
def test_gammaln():
    assert_warns_message(DeprecationWarning, "The function gammaln"
                         " is deprecated in 0.18 and will be removed"
                         " in 0.20. Use scipy.special.gammaln instead.",
                         gammaln, 3)


@ignore_warnings(category=DeprecationWarning)
def test_log_normalize():
    v = np.array([0.1, 0.8, 0.01, 0.09])
    a = np.log(2 * v)
    result = assert_warns_message(DeprecationWarning, "The function "
                                  "log_normalize is deprecated in 0.18 and"
                                  " will be removed in 0.20.",
                                  log_normalize, a)
    assert np.allclose(v, result, rtol=0.01)


@ignore_warnings(category=DeprecationWarning)
def test_wishart_log_det():
    a = np.array([0.1, 0.8, 0.01, 0.09])
    b = np.array([0.2, 0.7, 0.05, 0.1])
    assert_warns_message(DeprecationWarning, "The function "
                         "wishart_log_det is deprecated in 0.18 and"
                         " will be removed in 0.20.",
                         wishart_log_det, a, b, 2, 4)


@ignore_warnings(category=DeprecationWarning)
def test_wishart_logz():
    assert_warns_message(DeprecationWarning, "The function "
                         "wishart_logz is deprecated in 0.18 and "
                         "will be removed in 0.20.", wishart_logz,
                         3, np.identity(3), 1, 3)


@ignore_warnings(category=DeprecationWarning)
def test_DPGMM_deprecation():
    assert_warns_message(
      DeprecationWarning, "The `DPGMM` class is not working correctly and "
      "it's better to use `sklearn.mixture.BayesianGaussianMixture` class "
      "with parameter `weight_concentration_prior_type='dirichlet_process'` "
      "instead. DPGMM is deprecated in 0.18 and will be removed in 0.20.",
      DPGMM)


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


def test_VBGMM_deprecation():
    assert_warns_message(
        DeprecationWarning, "The `VBGMM` class is not working correctly and "
        "it's better to use `sklearn.mixture.BayesianGaussianMixture` class "
        "with parameter `weight_concentration_prior_type="
        "'dirichlet_distribution'` instead. VBGMM is deprecated "
        "in 0.18 and will be removed in 0.20.", VBGMM)


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


def test_vbgmm_no_modify_alpha():
    alpha = 2.
    n_components = 3
    X, y = make_blobs(random_state=1)
    vbgmm = VBGMM(n_components=n_components, alpha=alpha, n_iter=1)
    assert_equal(vbgmm.alpha, alpha)
    assert_equal(vbgmm.fit(X).alpha_, float(alpha) / n_components)
