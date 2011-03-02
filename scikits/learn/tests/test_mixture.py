import itertools
import unittest

import nose
from numpy.testing import assert_array_equal, assert_array_almost_equal, \
     assert_raises
import numpy as np
from scipy import stats

from scikits.learn import mixture

np.random.seed(0)


def _generate_random_spd_matrix(ndim):
    """Return a random symmetric, positive-definite matrix."""
    A = np.random.rand(ndim, ndim)
    U, s, V = np.linalg.svd(np.dot(A.T, A))
    randspd = np.dot(np.dot(U, 1.0 + np.diag(np.random.rand(ndim))), V)
    return randspd


def test_logsum_1D():
    A = np.random.rand(2) + 1.0
    for axis in range(1):
        Asum = mixture.logsum(A, axis)
        assert_array_almost_equal(np.exp(Asum), np.sum(np.exp(A), axis))


def test_logsum_3D():
    """
    Test also on a 3D matrix
    """
    A = np.random.rand(2, 2, 2) + 1.0
    for axis in range(3):
        Asum = mixture.logsum(A, axis)
        assert_array_almost_equal(np.exp(Asum), np.sum(np.exp(A), axis))


def test_normalize_1D():
    A = np.random.rand(2) + 1.0
    for axis in range(1):
        Anorm = mixture.normalize(A, axis)
        assert np.all(np.allclose(Anorm.sum(axis), 1.0))


def test_normalize_3D():
    A = np.random.rand(2, 2, 2) + 1.0
    for axis in range(3):
        Anorm = mixture.normalize(A, axis)
        assert np.all(np.allclose(Anorm.sum(axis), 1.0))


def test_sample_gaussian():
    """
    Test sample generation from mixture.sample_gaussian where covariance
    is diagonal, spherical and full
    """

    n_features, n_samples = 2, 300
    axis = 1
    mu = np.random.randint(10) * np.random.rand(n_features)
    cv = (np.random.rand(n_features) + 1.0) ** 2

    samples = mixture.sample_gaussian(
        mu, cv, cvtype='diag', n_samples=n_samples)

    assert np.allclose(samples.mean(axis), mu, atol=0.3)
    assert np.allclose(samples.var(axis),  cv, atol=0.5)

    # the same for spherical covariances
    cv = (np.random.rand() + 1.0) ** 2
    samples = mixture.sample_gaussian(
        mu, cv, cvtype='spherical', n_samples=n_samples)

    assert np.allclose(samples.mean(axis), mu, atol=0.3)
    assert np.allclose(
        samples.var(axis), np.repeat(cv, n_features), atol=0.5)

    # and for full covariances
    A = np.random.randn(n_features, n_features)
    cv = np.dot(A.T, A) + np.eye(n_features)
    samples = mixture.sample_gaussian(
        mu, cv, cvtype='full', n_samples=n_samples)
    assert np.allclose(samples.mean(axis), mu, atol=0.3)
    assert np.allclose(np.cov(samples), cv, atol=1.)


def _naive_lmvnpdf_diag(obs, mu, cv):
    # slow and naive implementation of lmvnpdf
    ref = np.empty((len(obs), len(mu)))
    stds = np.sqrt(cv)
    for i, (m, std) in enumerate(itertools.izip(mu, stds)):
        ref[:, i] = np.log(stats.norm.pdf(obs, m, std)).sum(axis=1)
    return ref


def test_lmvnpdf_diag():
    """
    test a slow and naive implementation of lmvnpdf and
    compare it to the vectorized version (mixture.lmvnpdf) to test
    for correctness
    """
    n_features, n_states, n_obs = 2, 3, 10
    mu = np.random.randint(10) * np.random.rand(n_states, n_features)
    cv = (np.random.rand(n_states, n_features) + 1.0) ** 2
    obs = np.random.randint(10) * np.random.rand(n_obs, n_features)

    ref = _naive_lmvnpdf_diag(obs, mu, cv)
    lpr = mixture.lmvnpdf(obs, mu, cv, 'diag')
    assert_array_almost_equal(lpr, ref)


def test_lmvnpdf_spherical():
    n_features, n_states, n_obs = 2, 3, 10

    mu = np.random.randint(10) * np.random.rand(n_states, n_features)
    spherecv = np.random.rand(n_states, 1) ** 2 + 1
    obs = np.random.randint(10) * np.random.rand(n_obs, n_features)

    cv = np.tile(spherecv, (n_features, 1))
    reference = _naive_lmvnpdf_diag(obs, mu, cv)
    lpr = mixture.lmvnpdf(obs, mu, spherecv, 'spherical')
    assert_array_almost_equal(lpr, reference)


def test_lmvnpdf_full():
    n_features, n_states, n_obs = 2, 3, 10

    mu = np.random.randint(10) * np.random.rand(n_states, n_features)
    cv = (np.random.rand(n_states, n_features) + 1.0) ** 2
    obs = np.random.randint(10) * np.random.rand(n_obs, n_features)

    fullcv = np.array([np.diag(x) for x in cv])

    reference = _naive_lmvnpdf_diag(obs, mu, cv)
    lpr = mixture.lmvnpdf(obs, mu, fullcv, 'full')
    assert_array_almost_equal(lpr, reference)


def test_GMM_attributes():
    n_states, n_features = 10, 4
    cvtype = 'diag'
    g = mixture.GMM(n_states, cvtype)
    weights = np.random.rand(n_states)
    weights = weights / weights.sum()
    means = np.random.randint(-20, 20, (n_states, n_features))

    assert g.n_states == n_states
    assert g.cvtype == cvtype

    g.weights = weights
    assert_array_almost_equal(g.weights, weights)
    assert_raises(ValueError, g.__setattr__, 'weights',
                      2 * weights)
    assert_raises(ValueError, g.__setattr__, 'weights', [])
    assert_raises(ValueError, g.__setattr__, 'weights',
                      np.zeros((n_states - 2, n_features)))

    g.means = means
    assert_array_almost_equal(g.means, means)
    assert_raises(ValueError, g.__setattr__, 'means', [])
    assert_raises(ValueError, g.__setattr__, 'means',
                      np.zeros((n_states - 2, n_features)))

    covars = (0.1 + 2 * np.random.rand(n_states, n_features)) ** 2
    g._covars = covars
    assert_array_almost_equal(g._covars, covars)
    assert_raises(ValueError, g.__setattr__, 'covars', [])
    assert_raises(ValueError, g.__setattr__, 'covars',
                      np.zeros((n_states - 2, n_features)))

    assert_raises(ValueError, mixture.GMM, n_states=20, cvtype='badcvtype')


class GMMTester():
    n_states = 10
    n_features = 4
    weights = np.random.rand(n_states)
    weights = weights / weights.sum()
    means = np.random.randint(-20, 20, (n_states, n_features))
    I = np.eye(n_features)
    covars = {'spherical': (0.1 + 2 * np.random.rand(n_states)) ** 2,
              'tied': _generate_random_spd_matrix(n_features) + 5 * I,
              'diag': (0.1 + 2 * np.random.rand(n_states, n_features)) ** 2,
              'full': np.array([_generate_random_spd_matrix(n_features) + 5 * I
                                for x in xrange(n_states)])}

    def test_eval(self):
        g = mixture.GMM(self.n_states, self.cvtype)
        # Make sure the means are far apart so posteriors.argmax()
        # picks the actual component used to generate the observations.
        g.means = 20 * self.means
        g._covars = self.covars[self.cvtype]
        g.weights = self.weights

        gaussidx = np.repeat(range(self.n_states), 5)
        nobs = len(gaussidx)
        obs = np.random.randn(nobs, self.n_features) + g.means[gaussidx]

        ll, posteriors = g.eval(obs)

        self.assertEqual(len(ll), nobs)
        self.assertEqual(posteriors.shape, (nobs, self.n_states))
        assert_array_almost_equal(posteriors.sum(axis=1), np.ones(nobs))
        assert_array_equal(posteriors.argmax(axis=1), gaussidx)

    def test_rvs(self, n=100):
        g = mixture.GMM(self.n_states, self.cvtype)
        # Make sure the means are far apart so posteriors.argmax()
        # picks the actual component used to generate the observations.
        g.means = 20 * self.means
        g._covars = np.maximum(self.covars[self.cvtype], 0.1)
        g.weights = self.weights

        samples = g.rvs(n)
        self.assertEquals(samples.shape, (n, self.n_features))

    def test_train(self, params='wmc'):
        g = mixture.GMM(self.n_states, self.cvtype)
        g.weights = self.weights
        g.means = self.means
        g._covars = 20 * self.covars[self.cvtype]

        # Create a training set by sampling from the predefined distribution.
        train_obs = g.rvs(n_samples=100)

        g.fit(train_obs, n_iter=0, init_params=params)

        # Do one training iteration at a time so we can keep track of
        # the log likelihood to make sure that it increases after each
        # iteration.
        trainll = []
        for iter in xrange(5):
            g.fit(train_obs, n_iter=1, params=params, init_params='',
                  min_covar=1e-1)
            trainll.append(g.score(train_obs).sum())
        # Note that the log likelihood will sometimes decrease by a
        # very small amount after it has more or less converged due to
        # the addition of min_covar to the covariance (to prevent
        # underflow).  This is why the threshold is set to -0.5
        # instead of 0.
        self.assertTrue(np.all(np.diff(trainll) > -0.5))


class TestGMMWithSphericalCovars(unittest.TestCase, GMMTester):
    cvtype = 'spherical'


class TestGMMWithDiagonalCovars(unittest.TestCase, GMMTester):
    cvtype = 'diag'


class TestGMMWithTiedCovars(unittest.TestCase, GMMTester):
    cvtype = 'tied'


class TestGMMWithFullCovars(unittest.TestCase, GMMTester):
    cvtype = 'full'


if __name__ == '__main__':
    nose.runmodule()
