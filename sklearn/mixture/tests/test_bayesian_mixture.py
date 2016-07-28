import numpy as np
from scipy import linalg
from scipy.special import gammaln, digamma

from sklearn import cluster
from sklearn.datasets import make_blobs
from sklearn.datasets.samples_generator import make_spd_matrix
from sklearn.utils.testing import assert_array_almost_equal
from sklearn.utils.testing import assert_raise_message
from sklearn.utils.testing import assert_almost_equal

from sklearn.mixture.bayesian_mixture import _log_dirichlet_norm
from sklearn.mixture.bayesian_mixture import _log_wishart_norm

from sklearn.mixture import BayesianGaussianMixture

from sklearn.mixture.gaussian_mixture import _estimate_gaussian_parameters


def test_log_dirichlet_norm():
    rng = np.random.RandomState(0)

    alpha = rng.rand(2)
    expected_norm = gammaln(np.sum(alpha)) - np.sum(gammaln(alpha))
    predected_norm = _log_dirichlet_norm(alpha)

    assert_almost_equal(expected_norm, predected_norm)


def test_log_wishart_norm():
    rng = np.random.RandomState(0)

    n_components, n_features = 5, 2
    nu = np.abs(rng.rand(n_components)) + 1.
    log_det_precisions_chol = n_features * np.log(range(2, 2 + n_components))

    expected_norm = np.empty(5)
    for k, (nu_k, log_det_k) in enumerate(zip(nu, log_det_precisions_chol)):
        expected_norm[k] = -(
            nu_k * (log_det_k + .5 * n_features * np.log(2.)) + np.sum(gammaln(
                .5 * (nu_k - np.arange(0, n_features)[:, np.newaxis])), 0))
    predected_norm = _log_wishart_norm(nu, log_det_precisions_chol, n_features)

    assert_almost_equal(expected_norm, predected_norm)


def test_bayesian_mixture_covariance_type():
    rng = np.random.RandomState(0)
    n_samples, n_features = 10, 2
    X = rng.rand(n_samples, n_features)

    covariance_type = 'bad_covariance_type'
    bgmm = BayesianGaussianMixture(covariance_type=covariance_type)
    assert_raise_message(ValueError,
                         "Invalid value for 'covariance_type': %s "
                         "'covariance_type' should be in "
                         "['spherical', 'tied', 'diag', 'full']"
                         % covariance_type,
                         bgmm.fit, X)


def test_bayesian_mixture_weights_prior_initialisation():
    rng = np.random.RandomState(0)
    n_samples, n_components, n_features = 10, 5, 2
    X = rng.rand(n_samples, n_features)

    # Check raise message for a bad value of alpha_init
    bad_alpha_init = 0.
    bgmm = BayesianGaussianMixture(alpha_init=bad_alpha_init)
    assert_raise_message(ValueError,
                         "The parameter 'alpha_init' should be "
                         "greater than 0., but got %.3f."
                         % bad_alpha_init,
                         bgmm.fit, X)

    # Check correct init for a given value of alpha_init
    alpha_init = rng.rand()
    bgmm = BayesianGaussianMixture(alpha_init=alpha_init).fit(X)
    assert_almost_equal(alpha_init, bgmm._alpha_prior)

    # Check correct init for the default value of alpha_init
    bgmm = BayesianGaussianMixture(n_components=n_components).fit(X)
    assert_almost_equal(1. / n_components, bgmm._alpha_prior)


def test_bayesian_mixture_means_prior_initialisation():
    rng = np.random.RandomState(0)
    n_samples, n_components, n_features = 10, 3, 2
    X = rng.rand(n_samples, n_features)

    # Check raise message for a bad value of beta_init
    bad_beta_init = 0.
    bgmm = BayesianGaussianMixture(beta_init=bad_beta_init)
    assert_raise_message(ValueError,
                         "The parameter 'beta_init' should be "
                         "greater than 0., but got %.3f."
                         % bad_beta_init,
                         bgmm.fit, X)

    # Check correct init for a given value of beta_init
    beta_init = rng.rand()
    bgmm = BayesianGaussianMixture(beta_init=beta_init).fit(X)
    assert_almost_equal(beta_init, bgmm._beta_prior)

    # Check correct init for the default value of beta_init
    bgmm = BayesianGaussianMixture().fit(X)
    assert_almost_equal(1., bgmm._beta_prior)

    # Check raise message for a bad shape of mean_init
    mean_init = rng.rand(n_features + 1)
    bgmm = BayesianGaussianMixture(n_components=n_components,
                                   mean_init=mean_init)
    assert_raise_message(ValueError,
                         "The parameter 'means' should have the shape of ",
                         bgmm.fit, X)

    # Check correct init for a given value of mean_init
    mean_init = rng.rand(n_features)
    bgmm = BayesianGaussianMixture(n_components=n_components,
                                   mean_init=mean_init).fit(X)
    assert_almost_equal(mean_init, bgmm._mean_prior)

    # Check correct init for the default value of bemean_initta
    bgmm = BayesianGaussianMixture(n_components=n_components).fit(X)
    assert_almost_equal(X.mean(axis=0), bgmm._mean_prior)


def test_bayesian_mixture_precisions_prior_initialisation():
    rng = np.random.RandomState(0)
    n_samples, n_features = 10, 2
    X = rng.rand(n_samples, n_features)

    # Check raise message for a bad value of nu_init
    bad_nu_init = n_features - 1.
    bgmm = BayesianGaussianMixture(nu_init=bad_nu_init)
    assert_raise_message(ValueError,
                         "The parameter 'nu_init' should be "
                         "greater than %d, but got %.3f."
                         % (n_features - 1, bad_nu_init),
                         bgmm.fit, X)

    # Check correct init for a given value of nu_init
    nu_init = rng.rand() + n_features - 1.
    bgmm = BayesianGaussianMixture(nu_init=nu_init).fit(X)
    assert_almost_equal(nu_init, bgmm._nu_prior)

    # Check correct init for the default value of nu_init
    nu_init_default = n_features
    bgmm = BayesianGaussianMixture(nu_init=nu_init_default).fit(X)
    assert_almost_equal(nu_init_default, bgmm._nu_prior)

    # Check correct init for a given value of covariance_init
    covariance_init = {
        'full': np.cov(X.T, bias=1),
        'tied': np.cov(X.T, bias=1),
        'diag': np.diag(np.atleast_2d(np.cov(X.T, bias=1))),
        'spherical': rng.rand()}

    bgmm = BayesianGaussianMixture()
    for cov_type in ['full', 'tied', 'diag', 'spherical']:
        print(cov_type)
        bgmm.covariance_type = cov_type
        bgmm.covariance_init = covariance_init[cov_type]
        bgmm.fit(X)
        assert_almost_equal(covariance_init[cov_type],
                            bgmm._covariance_prior)

    # Check raise message for a bad spherical value of covariance_init
    bad_covariance_init = -1.
    bgmm = BayesianGaussianMixture(covariance_type='spherical',
                                   covariance_init=bad_covariance_init)
    assert_raise_message(ValueError,
                         "The parameter 'spherical covariance_init' "
                         "should be greater than 0., but got %.3f."
                         % bad_covariance_init,
                         bgmm.fit, X)

    # Check correct init for the default value of covariance_init
    covariance_init_default = {
        'full': np.eye(X.shape[1]),
        'tied': np.eye(X.shape[1]),
        'diag': np.diag(np.atleast_2d(np.cov(X.T))),
        'spherical': np.var(X, axis=0).mean()}

    bgmm = BayesianGaussianMixture()
    for cov_type in ['full', 'tied', 'diag', 'spherical']:
        bgmm.covariance_type = cov_type
        bgmm.fit(X)
        assert_almost_equal(covariance_init_default[cov_type],
                            bgmm._covariance_prior)


def test_bayesian_mixture_check_is_fitted():
    rng = np.random.RandomState(0)
    n_samples, n_features = 10, 2

    # Check raise message
    bgmm = BayesianGaussianMixture()
    X = rng.rand(n_samples, n_features)
    assert_raise_message(ValueError,
                         'This BayesianGaussianMixture instance is not '
                         'fitted yet.', bgmm.score, X)


def test_bayesian_mixture_weights():
    rng = np.random.RandomState(0)
    n_samples, n_features = 10, 2

    X = rng.rand(n_samples, n_features)
    bgmm = BayesianGaussianMixture().fit(X)

    # Check the weights values
    expected_weights = bgmm.alpha_ / np.sum(bgmm.alpha_)
    predected_weights = bgmm.weights_

    assert_almost_equal(expected_weights, predected_weights)

    # Check the weights sum = 1
    assert_almost_equal(np.sum(bgmm.weights_), 1.0)

# ADD TEST FOR CHECK THE VALUES OF THE COVARIANCES + CHECK LOWER_BOUND IS ALWAYS POSITIVE
