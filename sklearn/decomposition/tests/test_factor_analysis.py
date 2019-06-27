# Author: Christian Osendorfer <osendorf@gmail.com>
#         Alexandre Gramfort <alexandre.gramfort@inria.fr>
# License: BSD3

import numpy as np

from sklearn.utils.testing import assert_warns
from sklearn.utils.testing import assert_equal
from sklearn.utils.testing import assert_greater
from sklearn.utils.testing import assert_less
from sklearn.utils.testing import assert_raises
from sklearn.utils.testing import assert_almost_equal
from sklearn.utils.testing import assert_array_almost_equal
from sklearn.exceptions import ConvergenceWarning
from sklearn.decomposition import FactorAnalysis
from sklearn.utils.testing import ignore_warnings


# Ignore warnings from switching to more power iterations in randomized_svd
@ignore_warnings
def test_factor_analysis():
    # Test FactorAnalysis ability to recover the data covariance structure
    rng = np.random.RandomState(0)
    n_samples, n_features, n_components = 20, 5, 3

    # Some random settings for the generative model
    W = rng.randn(n_components, n_features)
    # latent variable of dim 3, 20 of it
    h = rng.randn(n_samples, n_components)
    # using gamma to model different noise variance
    # per component
    noise = rng.gamma(1, size=n_features) * rng.randn(n_samples, n_features)

    # generate observations
    # wlog, mean is 0
    X = np.dot(h, W) + noise

    assert_raises(ValueError, FactorAnalysis, svd_method='foo')
    fa_fail = FactorAnalysis()
    fa_fail.svd_method = 'foo'
    assert_raises(ValueError, fa_fail.fit, X)
    fas = []
    for method in ['randomized', 'lapack']:
        fa = FactorAnalysis(n_components=n_components, svd_method=method)
        fa.fit(X)
        fas.append(fa)

        X_t = fa.transform(X)
        assert_equal(X_t.shape, (n_samples, n_components))

        assert_almost_equal(fa.loglike_[-1], fa.score_samples(X).sum())
        assert_almost_equal(fa.score_samples(X).mean(), fa.score(X))

        diff = np.all(np.diff(fa.loglike_))
        assert_greater(diff, 0., 'Log likelihood dif not increase')

        # Sample Covariance
        scov = np.cov(X, rowvar=0., bias=1.)

        # Model Covariance
        mcov = fa.get_covariance()
        diff = np.sum(np.abs(scov - mcov)) / W.size
        assert_less(diff, 0.1, "Mean absolute difference is %f" % diff)
        fa = FactorAnalysis(n_components=n_components,
                            noise_variance_init=np.ones(n_features))
        assert_raises(ValueError, fa.fit, X[:, :2])

    f = lambda x, y: np.abs(getattr(x, y))  # sign will not be equal
    fa1, fa2 = fas
    for attr in ['loglike_', 'components_', 'noise_variance_']:
        assert_almost_equal(f(fa1, attr), f(fa2, attr))

    fa1.max_iter = 1
    fa1.verbose = True
    assert_warns(ConvergenceWarning, fa1.fit, X)

    # Test get_covariance and get_precision with n_components == n_features
    # with n_components < n_features and with n_components == 0
    for n_components in [0, 2, X.shape[1]]:
        fa.n_components = n_components
        fa.fit(X)
        cov = fa.get_covariance()
        precision = fa.get_precision()
        assert_array_almost_equal(np.dot(cov, precision),
                                  np.eye(X.shape[1]), 12)
