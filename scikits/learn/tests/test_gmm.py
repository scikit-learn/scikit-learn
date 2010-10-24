import itertools
import unittest

from numpy.testing import assert_array_equal, assert_array_almost_equal
import numpy as np
from scipy import stats

from .. import gmm

np.random.seed(0)


def _generate_random_spd_matrix(ndim):
    """Return a random symmetric, positive-definite matrix."""
    A = np.random.rand(ndim, ndim)
    U, s, V = np.linalg.svd(np.dot(A.T, A))
    randspd = np.dot(np.dot(U, 1.0 + np.diag(np.random.rand(ndim))), V)
    return randspd


def test_simple1():
    X = [[0, 0], [.1, .1], [1, 1]]
    clf = gmm.GMM(2, n_dim=2)
    clf.fit(X)


class TestLogsum(unittest.TestCase):

    def test_logsum_1D(self):
        A = np.random.rand(10) + 1.0
        Asum = gmm.logsum(A)
        self.assertAlmostEqual(np.exp(Asum), np.sum(np.exp(A)))

    def test_logsum_2D(self):
        A = np.random.rand(10, 4) + 1.0
        Asum = gmm.logsum(A)
        self.assertAlmostEqual(np.exp(Asum), np.sum(np.exp(A)))

    def test_logsum_with_axis_1D(self):
        A = np.random.rand(10) + 1.0
        for axis in range(1):
            Asum = gmm.logsum(A, axis)
            assert_array_almost_equal(np.exp(Asum), np.sum(np.exp(A), axis))

    def test_logsum_with_axis_2D(self):
        A = np.random.rand(10, 4) + 1.0
        for axis in range(2):
            Asum = gmm.logsum(A, axis)
            assert_array_almost_equal(np.exp(Asum), np.sum(np.exp(A), axis))

    def test_logsum_with_axis_3D(self):
        A = np.random.rand(10, 4, 5) + 1.0
        for axis in range(3):
            Asum = gmm.logsum(A, axis)
            assert_array_almost_equal(np.exp(Asum), np.sum(np.exp(A), axis))


class TestNormalize(unittest.TestCase):

    def test_normalize_1D(self):
        A = np.random.rand(10) + 1.0
        Anorm = gmm.normalize(A)
        self.assertAlmostEqual(Anorm.sum(), 1.0)

    def test_normalize_2D(self):
        A = np.random.rand(10, 4) + 1.0
        Anorm = gmm.normalize(A)
        self.assertAlmostEqual(Anorm.sum(), 1.0)

    def test_normalize_with_axis_1D(self):
        A = np.random.rand(10) + 1.0
        for axis in range(1):
            Anorm = gmm.normalize(A, axis)
            self.assertTrue(np.all(np.allclose(Anorm.sum(axis), 1.0)))

    def test_normalize_with_axis_2D(self):
        A = np.random.rand(10, 4) + 1.0
        for axis in range(2):
            Anorm = gmm.normalize(A, axis)
            self.assertTrue(np.all(np.allclose(Anorm.sum(axis), 1.0)))

    def test_normalize_with_axis_3D(self):
        A = np.random.rand(10, 4, 5) + 1.0
        for axis in range(3):
            Anorm = gmm.normalize(A, axis)
            self.assertTrue(np.all(np.allclose(Anorm.sum(axis), 1.0)))


class TestSampleGaussian(unittest.TestCase):

    def _test_sample_gaussian_diag(self, n_dim, n=10000):
        mu = np.random.randint(10) * np.random.rand(n_dim)
        cv = (np.random.rand(n_dim) + 1.0) ** 2

        samples = gmm.sample_gaussian(mu, cv, cvtype='diag', n=n)

        if n_dim > 1:
            axis = 1
        else:
            axis = None
        assert_array_almost_equal(samples.mean(axis), mu, decimal=1)
        assert_array_almost_equal(samples.var(axis), cv, decimal=1)

    def test_sample_gaussian_diag_1D(self):
        self._test_sample_gaussian_diag(1)

    def test_sample_gaussian_diag_2D(self):
        self._test_sample_gaussian_diag(2)

    def test_sample_gaussian_diag_5D(self):
        self._test_sample_gaussian_diag(5)

    def _test_sample_gaussian_spherical(self, n_dim, n=10000):
        mu = np.random.randint(10) * np.random.rand(n_dim)
        cv = (np.random.rand() + 1.0) ** 2

        samples = gmm.sample_gaussian(mu, cv, cvtype='spherical', n=n)

        if n_dim > 1:
            axis = 1
        else:
            axis = None
        assert_array_almost_equal(samples.mean(axis), mu, decimal=1)
        assert_array_almost_equal(samples.var(axis), np.repeat(cv, n_dim),
                                  decimal=1)

    def test_sample_gaussian_spherical_1D(self):
        self._test_sample_gaussian_spherical(1)

    def test_sample_gaussian_spherical_2D(self):
        self._test_sample_gaussian_spherical(2)

    def test_sample_gaussian_spherical_5D(self):
        self._test_sample_gaussian_spherical(5)

    def _test_sample_gaussian_full(self, n_dim, n=10000):
        mu = np.random.randint(10) * np.random.rand(n_dim)
        cv = _generate_random_spd_matrix(n_dim) 

        samples = gmm.sample_gaussian(mu, cv, cvtype='full', n=n)

        if n_dim > 1:
            axis = 1
        else:
            axis = None
        assert_array_almost_equal(samples.mean(axis), mu, decimal=1)
        assert_array_almost_equal(np.cov(samples), cv, decimal=1)

    def test_sample_gaussian_full_1D(self):
        self._test_sample_gaussian_full(1)

    def test_sample_gaussian_full_2D(self):
        self._test_sample_gaussian_full(2)

    def test_sample_gaussian_full_5D(self):
        self._test_sample_gaussian_full(5)


class TestLmvnpdf(unittest.TestCase):

    def _slow_lmvnpdfdiag(self, obs, means, covars):
        lpr = np.empty((len(obs), len(means)))
        stds = np.sqrt(covars)
        for c, (mu, std) in enumerate(itertools.izip(means, stds)):
            for o, currobs in enumerate(obs):
                lpr[o,c] = np.log(stats.norm.pdf(currobs, mu, std)).sum()
        return lpr

    def _test_lmvnpdfdiag(self, n_dim, n_states, nobs=100):
        # test the slow and naive implementation of lmvnpdf and
        # compare it to the vectorized version (gmm.lmvnpdf) to test
        # for correctness
        mu = np.random.randint(10) * np.random.rand(n_states, n_dim)
        cv = (np.random.rand(n_states, n_dim) + 1.0) ** 2
        obs = np.random.randint(10) * np.random.rand(nobs, n_dim)

        reference = self._slow_lmvnpdfdiag(obs, mu, cv)
        lpr = gmm.lmvnpdf(obs, mu, cv, 'diag')
        assert_array_almost_equal(lpr, reference)

    def test_lmvnpdfdiag_univariate_single_gaussian(self):
        self._test_lmvnpdfdiag(1, 1)

    def test_lmvnpdfdiag_univariate(self):
        self._test_lmvnpdfdiag(1, 10)

    def test_lmvnpdfdiag_single_gaussian(self):
        self._test_lmvnpdfdiag(5, 1)

    def test_lmvnpdfdiag(self):
        self._test_lmvnpdfdiag(5, 10)

    def _test_lmvnpdfspherical(self, n_dim, n_states, nobs=100):
        mu = np.random.randint(10) * np.random.rand(n_states, n_dim)
        spherecv = np.random.rand(n_states, 1) ** 2 + 1
        obs = np.random.randint(10) * np.random.rand(nobs, n_dim)

        cv = np.tile(spherecv, (n_dim, 1))
        reference = self._slow_lmvnpdfdiag(obs, mu, cv)
        lpr = gmm.lmvnpdf(obs, mu, spherecv, 'spherical')
        assert_array_almost_equal(lpr, reference)

    def test_lmvnpdfspherical_univariate_single_gaussian(self):
        self._test_lmvnpdfspherical(1, 1)

    def test_lmvnpdfspherical_univariate(self):
        self._test_lmvnpdfspherical(1, 10)

    def test_lmvnpdfspherical_single_gaussian(self):
        self._test_lmvnpdfspherical(5, 1)

    def test_lmvnpdfspherical(self):
        self._test_lmvnpdfspherical(5, 10)

    def _test_lmvnpdffull_with_diagonal_covariance(self, n_dim, n_states,
                                                   nobs=100):
        mu = np.random.randint(10) * np.random.rand(n_states, n_dim)
        cv = (np.random.rand(n_states, n_dim) + 1.0) ** 2
        obs = np.random.randint(10) * np.random.rand(nobs, n_dim)

        fullcv = np.array([np.diag(x) for x in cv])

        reference = self._slow_lmvnpdfdiag(obs, mu, cv)
        lpr = gmm.lmvnpdf(obs, mu, fullcv, 'full')
        assert_array_almost_equal(lpr, reference)

    def test_lmvnpdffull_with_diagonal_covariance_univariate_single_gaussian(
        self):
        self._test_lmvnpdffull_with_diagonal_covariance(1, 1)

    def test_lmvnpdffull_with_diagonal_covariance_univariate(self):
        self._test_lmvnpdffull_with_diagonal_covariance(1, 10)

    def test_lmvnpdffull_with_diagonal_covariance_single_gaussian(self):
        self._test_lmvnpdffull_with_diagonal_covariance(5, 1)

    def test_lmvnpdffull_with_diagonal_covariance(self):
        self._test_lmvnpdffull_with_diagonal_covariance(5, 10)

    def _test_lmvnpdftied_with_diagonal_covariance(self, n_dim, n_states,
                                                   nobs=100):
        mu = np.random.randint(10) * np.random.rand(n_states, n_dim)
        tiedcv = (np.random.rand(n_dim) + 1.0) ** 2
        obs = np.random.randint(10) * np.random.rand(nobs, n_dim)

        cv = np.tile(tiedcv, (n_states, 1))

        reference = self._slow_lmvnpdfdiag(obs, mu, cv)
        lpr = gmm.lmvnpdf(obs, mu, np.diag(tiedcv), 'tied')
        assert_array_almost_equal(lpr, reference)

    def test_lmvnpdftied_with_diagonal_covariance_univariate_single_gaussian(
        self):
        self._test_lmvnpdftied_with_diagonal_covariance(1, 1)

    def test_lmvnpdftied_with_diagonal_covariance_univariate(self):
        self._test_lmvnpdftied_with_diagonal_covariance(1, 10)

    def test_lmvnpdftied_with_diagonal_covariance_single_gaussian(self):
        self._test_lmvnpdftied_with_diagonal_covariance(5, 1)

    def test_lmvnpdftied_with_diagonal_covariance(self):
        self._test_lmvnpdftied_with_diagonal_covariance(5, 10)

    def test_lmvnpdftied_consistent_with_lmvnpdffull(self):
        n_states = 4
        n_dim = 20
        nobs = 200

        mu = np.random.randint(10) * np.random.rand(n_states, n_dim)
        tiedcv = _generate_random_spd_matrix(n_dim)
        obs = np.random.randint(10) * np.random.rand(nobs, n_dim)

        cv = np.tile(tiedcv, (n_states, 1, 1))

        reference = gmm.lmvnpdf(obs, mu, cv, 'full')
        lpr = gmm.lmvnpdf(obs, mu, tiedcv, 'tied')
        assert_array_almost_equal(lpr, reference)


class GMMTester():
    n_states = 10
    n_dim = 4
    weights = np.random.rand(n_states)
    weights = weights / weights.sum()
    means = np.random.randint(-20, 20, (n_states, n_dim))
    I = np.eye(n_dim)
    covars = {'spherical': (0.1 + 2 * np.random.rand(n_states)) ** 2,
              'tied': _generate_random_spd_matrix(n_dim) + 5 * I,
              'diag': (0.1 + 2 * np.random.rand(n_states, n_dim)) ** 2,
              'full': np.array([_generate_random_spd_matrix(n_dim) + 5 * I
                                for x in xrange(n_states)])}

    def test_bad_cvtype(self):
        gmm.GMM(20, self.n_dim, self.cvtype)

        self.assertRaises(ValueError, gmm.GMM, 20, 1, 'badcvtype')

    def test_attributes(self):
        g = gmm.GMM(self.n_states, self.n_dim, self.cvtype)
        self.assertEquals(g.n_states, self.n_states)
        self.assertEquals(g.cvtype, self.cvtype)

        g.weights = self.weights
        assert_array_almost_equal(g.weights, self.weights)
        self.assertRaises(ValueError, g.__setattr__, 'weights',
                          2 * self.weights)
        self.assertRaises(ValueError, g.__setattr__, 'weights', [])
        self.assertRaises(ValueError, g.__setattr__, 'weights',
                          np.zeros((self.n_states - 2, self.n_dim)))

        g.means = self.means
        assert_array_almost_equal(g.means, self.means)
        self.assertRaises(ValueError, g.__setattr__, 'means', [])
        self.assertRaises(ValueError, g.__setattr__, 'means',
                          np.zeros((self.n_states - 2, self.n_dim)))

        g._covars = self.covars[self.cvtype]
        assert_array_almost_equal(g._covars, self.covars[self.cvtype])
        self.assertRaises(ValueError, g.__setattr__, 'covars', [])
        self.assertRaises(ValueError, g.__setattr__, 'covars',
                          np.zeros((self.n_states - 2, self.n_dim)))

    def test_eval(self):
        g = gmm.GMM(self.n_states, self.n_dim, self.cvtype)
        # Make sure the means are far apart so posteriors.argmax()
        # picks the actual component used to generate the observations.
        g.means = 20 * self.means
        g._covars = self.covars[self.cvtype]

        gaussidx = np.repeat(range(self.n_states), 5)
        nobs = len(gaussidx)
        obs = np.random.randn(nobs, self.n_dim) + g.means[gaussidx]

        ll, posteriors = g.eval(obs)

        self.assertEqual(len(ll), nobs)
        self.assertEqual(posteriors.shape, (nobs, self.n_states))
        assert_array_almost_equal(posteriors.sum(axis=1), np.ones(nobs))
        assert_array_equal(posteriors.argmax(axis=1), gaussidx)

    def test_rvs(self, n=1000):
        g = gmm.GMM(self.n_states, self.n_dim, self.cvtype)
        # Make sure the means are far apart so posteriors.argmax()
        # picks the actual component used to generate the observations.
        g.means = 20 * self.means
        g._covars = np.maximum(self.covars[self.cvtype], 0.1)
        # g.weights = self.weights

        samples = g.rvs(n)
        self.assertEquals(samples.shape, (n, self.n_dim))

    def test_train(self, params='wmc'):
        g = gmm.GMM(self.n_states, self.n_dim, self.cvtype)
        g.weights = self.weights
        g.means = self.means
        g._covars = 20 * self.covars[self.cvtype]

        # Create a training set by sampling from the predefined distribution.
        train_obs = g.rvs(n=100)

        g.fit(train_obs, n_iter=0, init_params=params)

        # Do one training iteration at a time so we can keep track of
        # the log likelihood to make sure that it increases after each
        # iteration.
        trainll = []
        for iter in xrange(20):
            g.fit(train_obs, n_iter=1, params=params, init_params='',
                  min_covar=1e-1)
            trainll.append(g.score(train_obs).sum())
        # Note that the log likelihood will sometimes decrease by a
        # very small amount after it has more or less converged due to
        # the addition of min_covar to the covariance (to prevent
        # underflow).  This is why the threshold is set to -0.5
        # instead of 0.
        print trainll, np.diff(trainll)
        self.assertTrue(np.all(np.diff(trainll) > -0.5))


class TestGMMWithSphericalCovars(unittest.TestCase, GMMTester):
    cvtype = 'spherical'


class TestGMMWithDiagonalCovars(unittest.TestCase, GMMTester):
    cvtype = 'diag'


class TestGMMWithTiedCovars(unittest.TestCase, GMMTester):
    cvtype = 'tied'


class TestGMMWithFullCovars(unittest.TestCase, GMMTester):
    cvtype = 'full'
