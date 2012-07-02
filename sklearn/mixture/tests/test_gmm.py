import itertools
import unittest

import numpy as np
from numpy.testing import assert_array_equal, assert_array_almost_equal, \
     assert_raises
from scipy import stats
from sklearn import mixture
from sklearn.datasets.samples_generator import make_spd_matrix

rng = np.random.RandomState(0)


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
        mu, cv, covariance_type='diag', n_samples=n_samples)

    assert np.allclose(samples.mean(axis), mu, atol=1.3)
    assert np.allclose(samples.var(axis), cv, atol=1.5)

    # the same for spherical covariances
    cv = (rng.rand() + 1.0) ** 2
    samples = mixture.sample_gaussian(
        mu, cv, covariance_type='spherical', n_samples=n_samples)

    assert np.allclose(samples.mean(axis), mu, atol=1.5)
    assert np.allclose(
        samples.var(axis), np.repeat(cv, n_features), atol=1.5)

    # and for full covariances
    A = rng.randn(n_features, n_features)
    cv = np.dot(A.T, A) + np.eye(n_features)
    samples = mixture.sample_gaussian(
        mu, cv, covariance_type='full', n_samples=n_samples)
    assert np.allclose(samples.mean(axis), mu, atol=1.3)
    assert np.allclose(np.cov(samples), cv, atol=2.5)


def _naive_lmvnpdf_diag(X, mu, cv):
    # slow and naive implementation of lmvnpdf
    ref = np.empty((len(X), len(mu)))
    stds = np.sqrt(cv)
    for i, (m, std) in enumerate(itertools.izip(mu, stds)):
        ref[:, i] = np.log(stats.norm.pdf(X, m, std)).sum(axis=1)
    return ref


def test_lmvnpdf_diag():
    """
    test a slow and naive implementation of lmvnpdf and
    compare it to the vectorized version (mixture.lmvnpdf) to test
    for correctness
    """
    n_features, n_components, n_samples = 2, 3, 10
    mu = rng.randint(10) * rng.rand(n_components, n_features)
    cv = (rng.rand(n_components, n_features) + 1.0) ** 2
    X = rng.randint(10) * rng.rand(n_samples, n_features)

    ref = _naive_lmvnpdf_diag(X, mu, cv)
    lpr = mixture.log_multivariate_normal_density(X, mu, cv, 'diag')
    assert_array_almost_equal(lpr, ref)


def test_lmvnpdf_spherical():
    n_features, n_components, n_samples = 2, 3, 10

    mu = rng.randint(10) * rng.rand(n_components, n_features)
    spherecv = rng.rand(n_components, 1) ** 2 + 1
    X = rng.randint(10) * rng.rand(n_samples, n_features)

    cv = np.tile(spherecv, (n_features, 1))
    reference = _naive_lmvnpdf_diag(X, mu, cv)
    lpr = mixture.log_multivariate_normal_density(X, mu, spherecv,
                                                  'spherical')
    assert_array_almost_equal(lpr, reference)


def test_lmvnpdf_full():
    n_features, n_components, n_samples = 2, 3, 10

    mu = rng.randint(10) * rng.rand(n_components, n_features)
    cv = (rng.rand(n_components, n_features) + 1.0) ** 2
    X = rng.randint(10) * rng.rand(n_samples, n_features)

    fullcv = np.array([np.diag(x) for x in cv])

    reference = _naive_lmvnpdf_diag(X, mu, cv)
    lpr = mixture.log_multivariate_normal_density(X, mu, fullcv, 'full')
    assert_array_almost_equal(lpr, reference)


def test_GMM_attributes():
    n_components, n_features = 10, 4
    covariance_type = 'diag'
    g = mixture.GMM(n_components, covariance_type, random_state=rng)
    weights = rng.rand(n_components)
    weights = weights / weights.sum()
    means = rng.randint(-20, 20, (n_components, n_features))

    assert g.n_components == n_components
    assert g.covariance_type == covariance_type

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
        self.covars = {'spherical': (0.1 + 2 * \
                        rng.rand(self.n_components, self.n_features)) ** 2,
                  'tied': make_spd_matrix(self.n_features, random_state=0) +\
                        5 * self.I,
                  'diag': (0.1 + 2 * rng.rand(self.n_components,\
                        self.n_features)) ** 2,
                  'full': np.array([make_spd_matrix(self.n_features,\
                        random_state=0)
                      + 5 * self.I for x in range(self.n_components)])}

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

        gaussidx = np.repeat(range(self.n_components), 5)
        n_samples = len(gaussidx)
        X = rng.randn(n_samples, self.n_features) + g.means_[gaussidx]

        ll, responsibilities = g.eval(X)

        self.assertEqual(len(ll), n_samples)
        self.assertEqual(responsibilities.shape,
                         (n_samples, self.n_components))
        assert_array_almost_equal(responsibilities.sum(axis=1),
                                  np.ones(n_samples))
        assert_array_equal(responsibilities.argmax(axis=1), gaussidx)

    def test_sample(self, n=100):
        g = self.model(n_components=self.n_components,
                       covariance_type=self.covariance_type, random_state=rng)
        # Make sure the means are far apart so responsibilities.argmax()
        # picks the actual component used to generate the observations.
        g.means_ = 20 * self.means
        g.covars_ = np.maximum(self.covars[self.covariance_type], 0.1)
        g.weights_ = self.weights

        samples = g.sample(n)
        self.assertEquals(samples.shape, (n, self.n_features))

    def test_train(self, params='wmc'):
        g = mixture.GMM(n_components=self.n_components,
                        covariance_type=self.covariance_type)
        g.weights_ = self.weights
        g.means_ = self.means
        g.covars_ = 20 * self.covars[self.covariance_type]

        # Create a training set by sampling from the predefined distribution.
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
        for iter in xrange(5):
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
        delta_min = np.diff(trainll).min()
        self.assertTrue(
            delta_min > self.threshold,
            "The min nll increase is %f which is lower than the admissible"
            " threshold of %f, for model %s. The likelihoods are %s."
            % (delta_min, self.threshold, self.covariance_type, trainll))

    def test_train_degenerate(self, params='wmc'):
        """ Train on degenerate data with 0 in some dimensions
        """
        # Create a training set by sampling from the predefined distribution.
        X = rng.randn(100, self.n_features)
        X.T[1:] = 0
        g = self.model(n_components=2, covariance_type=self.covariance_type,
                       random_state=rng, min_covar=1e-3, n_iter=5,
                       init_params=params)
        g.fit(X)
        trainll = g.score(X)
        self.assertTrue(np.sum(np.abs(trainll / 100 / X.shape[1])) < 5)

    def test_train_1d(self, params='wmc'):
        """ Train on 1-D data
        """
        # Create a training set by sampling from the predefined distribution.
        X = rng.randn(100, 1)
        #X.T[1:] = 0
        g = self.model(n_components=2, covariance_type=self.covariance_type,
                       random_state=rng, min_covar=1e-7, n_iter=5,
                       init_params=params)
        g.fit(X)
        trainll = g.score(X)
        if isinstance(g, mixture.DPGMM):
            self.assertTrue(np.sum(np.abs(trainll / 100)) < 5)
        else:
            self.assertTrue(np.sum(np.abs(trainll / 100)) < 2)

    def score(self, g, X):
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


def test_multiple_init():
    """Test that multiple inits does not much worse than a single one"""
    X = rng.randn(30, 5)
    X[:10] += 2
    g = mixture.GMM(n_components=2, covariance_type='spherical',
                    random_state=rng, min_covar=1e-7, n_iter=5)
    train1 = g.fit(X).score(X).sum()
    g.n_init = 5
    train2 = g.fit(X).score(X).sum()
    assert train2 >= train1 - 1.e-2


def test_n_parameters():
    """Test that the right number of parameters is estimated"""
    n_samples, n_dim, n_components = 7, 5, 2
    X = rng.randn(n_samples, n_dim)
    n_params = {'spherical': 13, 'diag': 21, 'tied': 26, 'full': 41}
    for cv_type in ['full', 'tied', 'diag', 'spherical']:
        g = mixture.GMM(n_components=n_components, covariance_type=cv_type,
                        random_state=rng, min_covar=1e-7, n_iter=1)
        g.fit(X)
        assert g._n_parameters() == n_params[cv_type]


def test_aic():
    """ Test the aic and bic criteria"""
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
        assert np.abs(g.aic(X) - aic) / n_samples < bound
        assert np.abs(g.bic(X) - bic) / n_samples < bound


if __name__ == '__main__':
    import nose
    nose.runmodule()
