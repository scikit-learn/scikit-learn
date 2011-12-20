import itertools
import unittest

import numpy as np
from numpy.testing import assert_array_equal, assert_array_almost_equal, \
     assert_raises
from scipy import stats

from sklearn import mixture
from sklearn.datasets.samples_generator import make_spd_matrix

rng = np.random.RandomState(0)


def test_normalize_1D():
    A = rng.rand(2) + 1.0
    for axis in range(1):
        Anorm = mixture.normalize(A, axis)
        assert np.all(np.allclose(Anorm.sum(axis), 1.0))


def test_normalize_3D():
    A = rng.rand(2, 2, 2) + 1.0
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
    mu = rng.randint(10) * rng.rand(n_features)
    cv = (rng.rand(n_features) + 1.0) ** 2

    samples = mixture.sample_gaussian(
        mu, cv, cvtype='diag', n_samples=n_samples)

    assert np.allclose(samples.mean(axis), mu, atol=1.3)
    assert np.allclose(samples.var(axis),  cv, atol=1.5)

    # the same for spherical covariances
    cv = (rng.rand() + 1.0) ** 2
    samples = mixture.sample_gaussian(
        mu, cv, cvtype='spherical', n_samples=n_samples)

    assert np.allclose(samples.mean(axis), mu, atol=1.5)
    assert np.allclose(
        samples.var(axis), np.repeat(cv, n_features), atol=1.5)

    # and for full covariances
    A = rng.randn(n_features, n_features)
    cv = np.dot(A.T, A) + np.eye(n_features)
    samples = mixture.sample_gaussian(
        mu, cv, cvtype='full', n_samples=n_samples)
    assert np.allclose(samples.mean(axis), mu, atol=1.3)
    assert np.allclose(np.cov(samples), cv, atol=2.5)


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
    n_features, n_components, n_obs = 2, 3, 10
    mu = rng.randint(10) * rng.rand(n_components, n_features)
    cv = (rng.rand(n_components, n_features) + 1.0) ** 2
    obs = rng.randint(10) * rng.rand(n_obs, n_features)

    ref = _naive_lmvnpdf_diag(obs, mu, cv)
    lpr = mixture.lmvnpdf(obs, mu, cv, 'diag')
    assert_array_almost_equal(lpr, ref)


def test_lmvnpdf_spherical():
    n_features, n_components, n_obs = 2, 3, 10

    mu = rng.randint(10) * rng.rand(n_components, n_features)
    spherecv = rng.rand(n_components, 1) ** 2 + 1
    obs = rng.randint(10) * rng.rand(n_obs, n_features)

    cv = np.tile(spherecv, (n_features, 1))
    reference = _naive_lmvnpdf_diag(obs, mu, cv)
    lpr = mixture.lmvnpdf(obs, mu, spherecv, 'spherical')
    assert_array_almost_equal(lpr, reference)


def test_lmvnpdf_full():
    n_features, n_components, n_obs = 2, 3, 10

    mu = rng.randint(10) * rng.rand(n_components, n_features)
    cv = (rng.rand(n_components, n_features) + 1.0) ** 2
    obs = rng.randint(10) * rng.rand(n_obs, n_features)

    fullcv = np.array([np.diag(x) for x in cv])

    reference = _naive_lmvnpdf_diag(obs, mu, cv)
    lpr = mixture.lmvnpdf(obs, mu, fullcv, 'full')
    assert_array_almost_equal(lpr, reference)


def test_GMM_attributes():
    n_components, n_features = 10, 4
    cvtype = 'diag'
    g = mixture.GMM(n_components, cvtype, random_state=rng)
    weights = rng.rand(n_components)
    weights = weights / weights.sum()
    means = rng.randint(-20, 20, (n_components, n_features))

    assert g.n_components == n_components
    assert g.cvtype == cvtype

    g.weights = weights
    assert_array_almost_equal(g.weights, weights)
    assert_raises(ValueError, g.__setattr__, 'weights',
                      2 * weights)
    assert_raises(ValueError, g.__setattr__, 'weights', [])
    assert_raises(ValueError, g.__setattr__, 'weights',
                      np.zeros((n_components - 2, n_features)))

    g.means = means
    assert_array_almost_equal(g.means, means)
    assert_raises(ValueError, g.__setattr__, 'means', [])
    assert_raises(ValueError, g.__setattr__, 'means',
                      np.zeros((n_components - 2, n_features)))

    covars = (0.1 + 2 * rng.rand(n_components, n_features)) ** 2
    g._covars = covars
    assert_array_almost_equal(g._covars, covars)
    assert_raises(ValueError, g.__setattr__, 'covars', [])
    assert_raises(ValueError, g.__setattr__, 'covars',
                      np.zeros((n_components - 2, n_features)))

    assert_raises(ValueError, mixture.GMM, n_components=20, cvtype='badcvtype')


class GMMTester():
    do_test_eval = True
    n_components = 10
    n_features = 4
    weights = rng.rand(n_components)
    weights = weights / weights.sum()
    means = rng.randint(-20, 20, (n_components, n_features))
    threshold = -0.5
    I = np.eye(n_features)
    covars = {'spherical': (0.1 + 2 * rng.rand(n_components)) ** 2,
              'tied': make_spd_matrix(n_features, random_state=0) + 5 * I,
              'diag': (0.1 + 2 * rng.rand(n_components, n_features)) ** 2,
              'full': np.array([make_spd_matrix(n_features, random_state=0)
                  + 5 * I for x in xrange(n_components)])}

    def test_eval(self):
        if not self.do_test_eval:
            return  # DPGMM does not support setting the means and
        # covariances before fitting There is no way of fixing this
        # due to the variational parameters being more expressive than
        # covariance matrices
        g = self.model(n_components=self.n_components, cvtype=self.cvtype,
                       random_state=rng)
        # Make sure the means are far apart so posteriors.argmax()
        # picks the actual component used to generate the observations.
        g.means = 20 * self.means
        g._covars = self.covars[self.cvtype]
        g.weights = self.weights

        gaussidx = np.repeat(range(self.n_components), 5)
        nobs = len(gaussidx)
        obs = rng.randn(nobs, self.n_features) + g.means[gaussidx]

        ll, posteriors = g.eval(obs)

        self.assertEqual(len(ll), nobs)
        self.assertEqual(posteriors.shape, (nobs, self.n_components))
        assert_array_almost_equal(posteriors.sum(axis=1), np.ones(nobs))
        assert_array_equal(posteriors.argmax(axis=1), gaussidx)

    def test_rvs(self, n=100):
        g = self.model(n_components=self.n_components, cvtype=self.cvtype,
                       random_state=rng)
        # Make sure the means are far apart so posteriors.argmax()
        # picks the actual component used to generate the observations.
        g.means = 20 * self.means
        g._covars = np.maximum(self.covars[self.cvtype], 0.1)
        g.weights = self.weights

        samples = g.rvs(n)
        self.assertEquals(samples.shape, (n, self.n_features))

    def test_train(self, params='wmc'):
        g = mixture.GMM(n_components=self.n_components, cvtype=self.cvtype)
        g.weights = self.weights
        g.means = self.means
        g._covars = 20 * self.covars[self.cvtype]

        # Create a training set by sampling from the predefined distribution.
        train_obs = g.rvs(n_samples=100)
        g = self.model(n_components=self.n_components, cvtype=self.cvtype,
                       random_state=rng, min_covar=1e-1)
        g.fit(train_obs, n_iter=1, init_params=params)

        # Do one training iteration at a time so we can keep track of
        # the log likelihood to make sure that it increases after each
        # iteration.
        trainll = []
        for iter in xrange(5):
            g.fit(train_obs, n_iter=1, params=params, init_params='')
            trainll.append(self.score(g, train_obs))
        g.fit(train_obs, n_iter=10, params=params, init_params='')  # finish
                                                                    # fitting

        # Note that the log likelihood will sometimes decrease by a
        # very small amount after it has more or less converged due to
        # the addition of min_covar to the covariance (to prevent
        # underflow).  This is why the threshold is set to -0.5
        # instead of 0.
        delta_min = np.diff(trainll).min()
        self.assertTrue(
            delta_min > self.threshold,
            "The min nll increase is %f which is lower than the admissible"
            " threshold of %f, for model %s. The likelihoods are %s."
                % (delta_min, self.threshold, self.cvtype, trainll))

    def score(self, g, train_obs):
        return g.score(train_obs).sum()


class TestGMMWithSphericalCovars(unittest.TestCase, GMMTester):
    cvtype = 'spherical'
    model = mixture.GMM


class TestGMMWithDiagonalCovars(unittest.TestCase, GMMTester):
    cvtype = 'diag'
    model = mixture.GMM


class TestGMMWithTiedCovars(unittest.TestCase, GMMTester):
    cvtype = 'tied'
    model = mixture.GMM


class TestGMMWithFullCovars(unittest.TestCase, GMMTester):
    cvtype = 'full'
    model = mixture.GMM


if __name__ == '__main__':
    import nose
    nose.runmodule()
