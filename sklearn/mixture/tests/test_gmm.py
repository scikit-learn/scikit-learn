# Important note for the deprecation cleaning of 0.20 :
# All the functions and classes of this file have been deprecated in 0.18.
# When you remove this file please remove the related files
# - 'sklearn/mixture/dpgmm.py'
# - 'sklearn/mixture/gmm.py'
# - 'sklearn/mixture/test_dpgmm.py'
import unittest
import copy
import sys

import numpy as np
from numpy.testing import (assert_array_equal, assert_array_almost_equal,
                           assert_raises)
from scipy import stats
from sklearn import mixture
from sklearn.datasets.samples_generator import make_spd_matrix
from sklearn.utils.testing import (assert_true, assert_greater,
                                   assert_raise_message, assert_warns_message,
                                   ignore_warnings)
from sklearn.metrics.cluster import adjusted_rand_score
from sklearn.externals.six.moves import cStringIO as StringIO


rng = np.random.RandomState(0)


def test_sample_gaussian():
    # Test sample generation from mixture.sample_gaussian where covariance
    # is diagonal, spherical and full

    n_features, n_samples = 2, 300
    axis = 1
    mu = rng.randint(10) * rng.rand(n_features)
    cv = (rng.rand(n_features) + 1.0) ** 2

    samples = mixture.gmm._sample_gaussian(
        mu, cv, covariance_type='diag', n_samples=n_samples)

    assert_true(np.allclose(samples.mean(axis), mu, atol=1.3))
    assert_true(np.allclose(samples.var(axis), cv, atol=1.5))

    # the same for spherical covariances
    cv = (rng.rand() + 1.0) ** 2
    samples = mixture.gmm._sample_gaussian(
        mu, cv, covariance_type='spherical', n_samples=n_samples)

    assert_true(np.allclose(samples.mean(axis), mu, atol=1.5))
    assert_true(np.allclose(
        samples.var(axis), np.repeat(cv, n_features), atol=1.5))

    # and for full covariances
    A = rng.randn(n_features, n_features)
    cv = np.dot(A.T, A) + np.eye(n_features)
    samples = mixture.gmm._sample_gaussian(
        mu, cv, covariance_type='full', n_samples=n_samples)
    assert_true(np.allclose(samples.mean(axis), mu, atol=1.3))
    assert_true(np.allclose(np.cov(samples), cv, atol=2.5))

    # Numerical stability check: in SciPy 0.12.0 at least, eigh may return
    # tiny negative values in its second return value.
    x = mixture.gmm._sample_gaussian(
        [0, 0], [[4, 3], [1, .1]], covariance_type='full', random_state=42)
    assert_true(np.isfinite(x).all())


def _naive_lmvnpdf_diag(X, mu, cv):
    # slow and naive implementation of lmvnpdf
    ref = np.empty((len(X), len(mu)))
    stds = np.sqrt(cv)
    for i, (m, std) in enumerate(zip(mu, stds)):
        ref[:, i] = np.log(stats.norm.pdf(X, m, std)).sum(axis=1)
    return ref


def test_lmvnpdf_diag():
    # test a slow and naive implementation of lmvnpdf and
    # compare it to the vectorized version (mixture.lmvnpdf) to test
    # for correctness
    n_features, n_components, n_samples = 2, 3, 10
    mu = rng.randint(10) * rng.rand(n_components, n_features)
    cv = (rng.rand(n_components, n_features) + 1.0) ** 2
    X = rng.randint(10) * rng.rand(n_samples, n_features)

    ref = _naive_lmvnpdf_diag(X, mu, cv)
    lpr = assert_warns_message(DeprecationWarning, "The function"
                             " log_multivariate_normal_density is "
                             "deprecated in 0.18 and will be removed in 0.20.",
                             mixture.log_multivariate_normal_density,
                             X, mu, cv, 'diag')
    assert_array_almost_equal(lpr, ref)


def test_lmvnpdf_spherical():
    n_features, n_components, n_samples = 2, 3, 10

    mu = rng.randint(10) * rng.rand(n_components, n_features)
    spherecv = rng.rand(n_components, 1) ** 2 + 1
    X = rng.randint(10) * rng.rand(n_samples, n_features)

    cv = np.tile(spherecv, (n_features, 1))
    reference = _naive_lmvnpdf_diag(X, mu, cv)
    lpr = assert_warns_message(DeprecationWarning, "The function"
                             " log_multivariate_normal_density is "
                             "deprecated in 0.18 and will be removed in 0.20.",
                             mixture.log_multivariate_normal_density,
                             X, mu, spherecv, 'spherical')
    assert_array_almost_equal(lpr, reference)

def test_lmvnpdf_full():
    n_features, n_components, n_samples = 2, 3, 10

    mu = rng.randint(10) * rng.rand(n_components, n_features)
    cv = (rng.rand(n_components, n_features) + 1.0) ** 2
    X = rng.randint(10) * rng.rand(n_samples, n_features)

    fullcv = np.array([np.diag(x) for x in cv])

    reference = _naive_lmvnpdf_diag(X, mu, cv)
    lpr = assert_warns_message(DeprecationWarning, "The function"
                             " log_multivariate_normal_density is "
                             "deprecated in 0.18 and will be removed in 0.20.",
                             mixture.log_multivariate_normal_density,
                             X, mu, fullcv, 'full')
    assert_array_almost_equal(lpr, reference)


def test_lvmpdf_full_cv_non_positive_definite():
    n_features, n_samples = 2, 10
    rng = np.random.RandomState(0)
    X = rng.randint(10) * rng.rand(n_samples, n_features)
    mu = np.mean(X, 0)
    cv = np.array([[[-1, 0], [0, 1]]])
    expected_message = "'covars' must be symmetric, positive-definite"
    assert_raise_message(ValueError, expected_message,
                         mixture.log_multivariate_normal_density,
                         X, mu, cv, 'full')


# This function tests the deprecated old GMM class
@ignore_warnings(category=DeprecationWarning)
def test_GMM_attributes():
    n_components, n_features = 10, 4
    covariance_type = 'diag'
    g = mixture.GMM(n_components, covariance_type, random_state=rng)
    weights = rng.rand(n_components)
    weights = weights / weights.sum()
    means = rng.randint(-20, 20, (n_components, n_features))

    assert_true(g.n_components == n_components)
    assert_true(g.covariance_type == covariance_type)

    g.weights_ = weights
    assert_array_almost_equal(g.weights_, weights)
    g.means_ = means
    assert_array_almost_equal(g.means_, means)

    covars = (0.1 + 2 * rng.rand(n_components, n_features)) ** 2
    g.covars_ = covars
    assert_array_almost_equal(g.covars_, covars)
    assert_raises(ValueError, g._set_covars, [])
    assert_raises(ValueError, g._set_covars,
                  np.zeros((n_components - 2, n_features)))

    assert_raises(ValueError, mixture.GMM, n_components=20,
                  covariance_type='badcovariance_type')


class GMMTester():
    do_test_eval = True

    def _setUp(self):
        self.n_components = 10
        self.n_features = 4
        self.weights = rng.rand(self.n_components)
        self.weights = self.weights / self.weights.sum()
        self.means = rng.randint(-20, 20, (self.n_components, self.n_features))
        self.threshold = -0.5
        self.I = np.eye(self.n_features)
        self.covars = {
            'spherical': (0.1 + 2 * rng.rand(self.n_components,
                                             self.n_features)) ** 2,
            'tied': (make_spd_matrix(self.n_features, random_state=0)
                     + 5 * self.I),
            'diag': (0.1 + 2 * rng.rand(self.n_components,
                                        self.n_features)) ** 2,
            'full': np.array([make_spd_matrix(self.n_features, random_state=0)
                              + 5 * self.I for x in range(self.n_components)])}

    # This function tests the deprecated old GMM class
    @ignore_warnings(category=DeprecationWarning)
    def test_eval(self):
        if not self.do_test_eval:
            return  # DPGMM does not support setting the means and
        # covariances before fitting There is no way of fixing this
        # due to the variational parameters being more expressive than
        # covariance matrices
        g = self.model(n_components=self.n_components,
                       covariance_type=self.covariance_type, random_state=rng)
        # Make sure the means are far apart so responsibilities.argmax()
        # picks the actual component used to generate the observations.
        g.means_ = 20 * self.means
        g.covars_ = self.covars[self.covariance_type]
        g.weights_ = self.weights

        gaussidx = np.repeat(np.arange(self.n_components), 5)
        n_samples = len(gaussidx)
        X = rng.randn(n_samples, self.n_features) + g.means_[gaussidx]

        with ignore_warnings(category=DeprecationWarning):
            ll, responsibilities = g.score_samples(X)

        self.assertEqual(len(ll), n_samples)
        self.assertEqual(responsibilities.shape,
                         (n_samples, self.n_components))
        assert_array_almost_equal(responsibilities.sum(axis=1),
                                  np.ones(n_samples))
        assert_array_equal(responsibilities.argmax(axis=1), gaussidx)

    # This function tests the deprecated old GMM class
    @ignore_warnings(category=DeprecationWarning)
    def test_sample(self, n=100):
        g = self.model(n_components=self.n_components,
                       covariance_type=self.covariance_type,
                       random_state=rng)
        # Make sure the means are far apart so responsibilities.argmax()
        # picks the actual component used to generate the observations.
        g.means_ = 20 * self.means
        g.covars_ = np.maximum(self.covars[self.covariance_type], 0.1)
        g.weights_ = self.weights

        with ignore_warnings(category=DeprecationWarning):
            samples = g.sample(n)
        self.assertEqual(samples.shape, (n, self.n_features))

    # This function tests the deprecated old GMM class
    @ignore_warnings(category=DeprecationWarning)
    def test_train(self, params='wmc'):
        g = mixture.GMM(n_components=self.n_components,
                        covariance_type=self.covariance_type)
        with ignore_warnings(category=DeprecationWarning):
            g.weights_ = self.weights
            g.means_ = self.means
            g.covars_ = 20 * self.covars[self.covariance_type]

        # Create a training set by sampling from the predefined distribution.
        with ignore_warnings(category=DeprecationWarning):
            X = g.sample(n_samples=100)
            g = self.model(n_components=self.n_components,
                           covariance_type=self.covariance_type,
                           random_state=rng, min_covar=1e-1,
                           n_iter=1, init_params=params)
            g.fit(X)

        # Do one training iteration at a time so we can keep track of
        # the log likelihood to make sure that it increases after each
        # iteration.
        trainll = []
        with ignore_warnings(category=DeprecationWarning):
            for _ in range(5):
                g.params = params
                g.init_params = ''
                g.fit(X)
                trainll.append(self.score(g, X))
            g.n_iter = 10
            g.init_params = ''
            g.params = params
            g.fit(X)  # finish fitting

        # Note that the log likelihood will sometimes decrease by a
        # very small amount after it has more or less converged due to
        # the addition of min_covar to the covariance (to prevent
        # underflow).  This is why the threshold is set to -0.5
        # instead of 0.
        with ignore_warnings(category=DeprecationWarning):
            delta_min = np.diff(trainll).min()
        self.assertTrue(
            delta_min > self.threshold,
            "The min nll increase is %f which is lower than the admissible"
            " threshold of %f, for model %s. The likelihoods are %s."
            % (delta_min, self.threshold, self.covariance_type, trainll))

    # This function tests the deprecated old GMM class
    @ignore_warnings(category=DeprecationWarning)
    def test_train_degenerate(self, params='wmc'):
        # Train on degenerate data with 0 in some dimensions
        # Create a training set by sampling from the predefined
        # distribution.
        X = rng.randn(100, self.n_features)
        X.T[1:] = 0
        g = self.model(n_components=2,
                       covariance_type=self.covariance_type,
                       random_state=rng, min_covar=1e-3, n_iter=5,
                       init_params=params)
        with ignore_warnings(category=DeprecationWarning):
            g.fit(X)
            trainll = g.score(X)
        self.assertTrue(np.sum(np.abs(trainll / 100 / X.shape[1])) < 5)

    # This function tests the deprecated old GMM class
    @ignore_warnings(category=DeprecationWarning)
    def test_train_1d(self, params='wmc'):
        # Train on 1-D data
        # Create a training set by sampling from the predefined
        # distribution.
        X = rng.randn(100, 1)
        # X.T[1:] = 0
        g = self.model(n_components=2,
                       covariance_type=self.covariance_type,
                       random_state=rng, min_covar=1e-7, n_iter=5,
                       init_params=params)
        with ignore_warnings(category=DeprecationWarning):
            g.fit(X)
            trainll = g.score(X)
            if isinstance(g, mixture.dpgmm._DPGMMBase):
                self.assertTrue(np.sum(np.abs(trainll / 100)) < 5)
            else:
                self.assertTrue(np.sum(np.abs(trainll / 100)) < 2)

    # This function tests the deprecated old GMM class
    @ignore_warnings(category=DeprecationWarning)
    def score(self, g, X):
        with ignore_warnings(category=DeprecationWarning):
            return g.score(X).sum()


class TestGMMWithSphericalCovars(unittest.TestCase, GMMTester):
    covariance_type = 'spherical'
    model = mixture.GMM
    setUp = GMMTester._setUp


class TestGMMWithDiagonalCovars(unittest.TestCase, GMMTester):
    covariance_type = 'diag'
    model = mixture.GMM
    setUp = GMMTester._setUp


class TestGMMWithTiedCovars(unittest.TestCase, GMMTester):
    covariance_type = 'tied'
    model = mixture.GMM
    setUp = GMMTester._setUp


class TestGMMWithFullCovars(unittest.TestCase, GMMTester):
    covariance_type = 'full'
    model = mixture.GMM
    setUp = GMMTester._setUp


# This function tests the deprecated old GMM class
@ignore_warnings(category=DeprecationWarning)
def test_multiple_init():
    # Test that multiple inits does not much worse than a single one
    X = rng.randn(30, 5)
    X[:10] += 2
    g = mixture.GMM(n_components=2, covariance_type='spherical',
                    random_state=rng, min_covar=1e-7, n_iter=5)
    with ignore_warnings(category=DeprecationWarning):
        train1 = g.fit(X).score(X).sum()
        g.n_init = 5
        train2 = g.fit(X).score(X).sum()
    assert_true(train2 >= train1 - 1.e-2)


# This function tests the deprecated old GMM class
@ignore_warnings(category=DeprecationWarning)
def test_n_parameters():
    n_samples, n_dim, n_components = 7, 5, 2
    X = rng.randn(n_samples, n_dim)
    n_params = {'spherical': 13, 'diag': 21, 'tied': 26, 'full': 41}
    for cv_type in ['full', 'tied', 'diag', 'spherical']:
        with ignore_warnings(category=DeprecationWarning):
            g = mixture.GMM(n_components=n_components, covariance_type=cv_type,
                            random_state=rng, min_covar=1e-7, n_iter=1)
            g.fit(X)
            assert_true(g._n_parameters() == n_params[cv_type])


# This function tests the deprecated old GMM class
@ignore_warnings(category=DeprecationWarning)
def test_1d_1component():
    # Test all of the covariance_types return the same BIC score for
    # 1-dimensional, 1 component fits.
    n_samples, n_dim, n_components = 100, 1, 1
    X = rng.randn(n_samples, n_dim)
    g_full = mixture.GMM(n_components=n_components, covariance_type='full',
                         random_state=rng, min_covar=1e-7, n_iter=1)
    with ignore_warnings(category=DeprecationWarning):
        g_full.fit(X)
        g_full_bic = g_full.bic(X)
        for cv_type in ['tied', 'diag', 'spherical']:
            g = mixture.GMM(n_components=n_components, covariance_type=cv_type,
                            random_state=rng, min_covar=1e-7, n_iter=1)
            g.fit(X)
            assert_array_almost_equal(g.bic(X), g_full_bic)


def assert_fit_predict_correct(model, X):
    model2 = copy.deepcopy(model)

    predictions_1 = model.fit(X).predict(X)
    predictions_2 = model2.fit_predict(X)

    assert adjusted_rand_score(predictions_1, predictions_2) == 1.0


# This function tests the deprecated old GMM class
@ignore_warnings(category=DeprecationWarning)
def test_fit_predict():
    """
    test that gmm.fit_predict is equivalent to gmm.fit + gmm.predict
    """
    lrng = np.random.RandomState(101)

    n_samples, n_dim, n_comps = 100, 2, 2
    mu = np.array([[8, 8]])
    component_0 = lrng.randn(n_samples, n_dim)
    component_1 = lrng.randn(n_samples, n_dim) + mu
    X = np.vstack((component_0, component_1))

    for m_constructor in (mixture.GMM, mixture.VBGMM, mixture.DPGMM):
        model = m_constructor(n_components=n_comps, covariance_type='full',
                              min_covar=1e-7, n_iter=5,
                              random_state=np.random.RandomState(0))
        assert_fit_predict_correct(model, X)

    model = mixture.GMM(n_components=n_comps, n_iter=0)
    z = model.fit_predict(X)
    assert np.all(z == 0), "Quick Initialization Failed!"


# This function tests the deprecated old GMM class
@ignore_warnings(category=DeprecationWarning)
def test_aic():
    # Test the aic and bic criteria
    n_samples, n_dim, n_components = 50, 3, 2
    X = rng.randn(n_samples, n_dim)
    SGH = 0.5 * (X.var() + np.log(2 * np.pi))  # standard gaussian entropy

    for cv_type in ['full', 'tied', 'diag', 'spherical']:
        g = mixture.GMM(n_components=n_components, covariance_type=cv_type,
                        random_state=rng, min_covar=1e-7)
        g.fit(X)
        aic = 2 * n_samples * SGH * n_dim + 2 * g._n_parameters()
        bic = (2 * n_samples * SGH * n_dim +
               np.log(n_samples) * g._n_parameters())
        bound = n_dim * 3. / np.sqrt(n_samples)
        assert_true(np.abs(g.aic(X) - aic) / n_samples < bound)
        assert_true(np.abs(g.bic(X) - bic) / n_samples < bound)


# This function tests the deprecated old GMM class
@ignore_warnings(category=DeprecationWarning)
def check_positive_definite_covars(covariance_type):
    r"""Test that covariance matrices do not become non positive definite

    Due to the accumulation of round-off errors, the computation of the
    covariance  matrices during the learning phase could lead to non-positive
    definite covariance matrices. Namely the use of the formula:

    .. math:: C = (\sum_i w_i  x_i x_i^T) - \mu \mu^T

    instead of:

    .. math:: C = \sum_i w_i (x_i - \mu)(x_i - \mu)^T

    while mathematically equivalent, was observed a ``LinAlgError`` exception,
    when computing a ``GMM`` with full covariance matrices and fixed mean.

    This function ensures that some later optimization will not introduce the
    problem again.
    """
    rng = np.random.RandomState(1)
    # we build a dataset with 2 2d component. The components are unbalanced
    # (respective weights 0.9 and 0.1)
    X = rng.randn(100, 2)
    X[-10:] += (3, 3)  # Shift the 10 last points

    gmm = mixture.GMM(2, params="wc", covariance_type=covariance_type,
                      min_covar=1e-3)

    # This is a non-regression test for issue #2640. The following call used
    # to trigger:
    # numpy.linalg.linalg.LinAlgError: 2-th leading minor not positive definite
    gmm.fit(X)

    if covariance_type == "diag" or covariance_type == "spherical":
        assert_greater(gmm.covars_.min(), 0)
    else:
        if covariance_type == "tied":
            covs = [gmm.covars_]
        else:
            covs = gmm.covars_

        for c in covs:
            assert_greater(np.linalg.det(c), 0)


def test_positive_definite_covars():
    # Check positive definiteness for all covariance types
    for covariance_type in ["full", "tied", "diag", "spherical"]:
        yield check_positive_definite_covars, covariance_type


# This function tests the deprecated old GMM class
@ignore_warnings(category=DeprecationWarning)
def test_verbose_first_level():
    # Create sample data
    X = rng.randn(30, 5)
    X[:10] += 2
    g = mixture.GMM(n_components=2, n_init=2, verbose=1)

    old_stdout = sys.stdout
    sys.stdout = StringIO()
    try:
        g.fit(X)
    finally:
        sys.stdout = old_stdout


# This function tests the deprecated old GMM class
@ignore_warnings(category=DeprecationWarning)
def test_verbose_second_level():
    # Create sample data
    X = rng.randn(30, 5)
    X[:10] += 2
    g = mixture.GMM(n_components=2, n_init=2, verbose=2)

    old_stdout = sys.stdout
    sys.stdout = StringIO()
    try:
        g.fit(X)
    finally:
        sys.stdout = old_stdout
