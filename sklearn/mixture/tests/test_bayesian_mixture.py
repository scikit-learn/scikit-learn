import numpy as np
from scipy import linalg
from scipy.special import gammaln, digamma

from sklearn.datasets.samples_generator import make_spd_matrix
from sklearn.utils.testing import assert_array_almost_equal
from sklearn.utils.testing import assert_raise_message
from sklearn.utils.testing import assert_almost_equal

from sklearn.mixture.bayesian_mixture import log_dirichlet_norm
from sklearn.mixture.bayesian_mixture import log_wishart_norm
from sklearn.mixture.bayesian_mixture import estimate_wishart_entropy
from sklearn.mixture.bayesian_mixture import gamma_entropy_spherical
from sklearn.mixture.bayesian_mixture import gamma_entropy_diag

from sklearn.mixture import BayesianGaussianMixture
from sklearn.datasets.samples_generator import make_spd_matrix

def test_log_dirichlet_norm():
    rng = np.random.RandomState(0)

    alpha = rng.rand(2)
    expected_norm = gammaln(np.sum(alpha)) - np.sum(gammaln(alpha))
    predected_norm = log_dirichlet_norm(alpha)

    assert_almost_equal(expected_norm, predected_norm)


def test_log_wishart_norm():
    rng = np.random.RandomState(0)

    n_features = 2
    nu = np.abs(rng.rand()) + 1.
    inv_W = linalg.inv(make_spd_matrix(n_features, rng))
    inv_W_chol = linalg.cholesky(inv_W, lower=True)

    expected_norm = (nu * np.sum(np.log(np.diag(inv_W_chol))) -
                     .5 * n_features * nu * np.log(2.) -
                     .25 * n_features * (n_features - 1) * np.log(np.pi) -
                     np.sum(gammaln(.5 * (nu + 1. -
                                          np.arange(1, n_features + 1.)))))
    predected_norm = log_wishart_norm(nu, inv_W_chol, n_features)

    assert_almost_equal(expected_norm, predected_norm)


def test_estimate_wishart_entropy():
    rng = np.random.RandomState(0)

    n_features = 2
    nu = np.abs(rng.rand()) + 1.
    inv_W = linalg.inv(make_spd_matrix(n_features, rng))
    inv_W_chol = linalg.cholesky(inv_W, lower=True)
    log_lambda = rng.rand()

    expected_entropy = (.5 * nu * n_features -
                        .5 * (nu - n_features - 1.) * log_lambda -
                        log_wishart_norm(nu, inv_W_chol, n_features))
    predected_entropy = estimate_wishart_entropy(nu, inv_W_chol, log_lambda,
                                                 n_features)

    assert_almost_equal(expected_entropy, predected_entropy)


def test_gamma_entropy_spherical():
    rng = np.random.RandomState(0)

    n_components = 5
    a = rng.rand(n_components)
    inv_b = rng.rand(n_components)

    expected_entropy = gammaln(a) - (a - 1.) * digamma(a) - np.log(inv_b) + a
    predected_entropy = gamma_entropy_spherical(a, inv_b)

    assert_almost_equal(expected_entropy, predected_entropy)


def test_gamma_entropy_diag():
    rng = np.random.RandomState(0)

    n_components, n_features = 5, 2
    a = rng.rand(n_components)
    inv_b = rng.rand(n_components, n_features)

    expected_entropy = ((gammaln(a) - (a - 1.) * digamma(a) + a) * len(inv_b) -
                        np.sum(np.log(inv_b)))
    predected_entropy = gamma_entropy_diag(a, inv_b)

    assert_almost_equal(expected_entropy, predected_entropy)


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

    # Check raise message for a bad value of alpha_prior_init
    bad_alpha_prior_init = 0.
    bgmm = BayesianGaussianMixture(alpha_prior_init=bad_alpha_prior_init)
    assert_raise_message(ValueError,
                         "The parameter 'alpha_prior_init' should be "
                         "greater than 0., but got %.3f."
                         % bad_alpha_prior_init,
                         bgmm.fit, X)

    # Check correct init for a given value of alpha_prior_init
    alpha_prior_init = rng.rand()
    bgmm = BayesianGaussianMixture(alpha_prior_init=alpha_prior_init).fit(X)
    assert_almost_equal(alpha_prior_init, bgmm._alpha_prior)

    # Check correct init for the default value of alpha_prior_init
    bgmm = BayesianGaussianMixture(n_components=n_components).fit(X)
    assert_almost_equal(1. / n_components, bgmm._alpha_prior)


def test_bayesian_mixture_means_prior_initialisation():
    rng = np.random.RandomState(0)
    n_samples, n_components, n_features = 10, 3, 2
    X = rng.rand(n_samples, n_features)

    # Check raise message for a bad value of beta_prior_init
    bad_beta_prior_init = 0.
    bgmm = BayesianGaussianMixture(beta_prior_init=bad_beta_prior_init)
    assert_raise_message(ValueError,
                         "The parameter 'beta_prior_init' should be "
                         "greater than 0., but got %.3f."
                         % bad_beta_prior_init,
                         bgmm.fit, X)

    # Check correct init for a given value of beta_prior_init
    beta_prior_init = rng.rand()
    bgmm = BayesianGaussianMixture(beta_prior_init=beta_prior_init).fit(X)
    assert_almost_equal(beta_prior_init, bgmm._beta_prior)

    # Check correct init for the default value of beta_prior_init
    bgmm = BayesianGaussianMixture().fit(X)
    assert_almost_equal(1., bgmm._beta_prior)

    # Check raise message for a bad shape of m_prior_init
    m_prior_init = rng.rand(n_features + 1)
    bgmm = BayesianGaussianMixture(n_components=n_components,
                                   m_prior_init=m_prior_init)
    assert_raise_message(ValueError,
                         "The parameter 'means' should have the shape of ",
                         bgmm.fit, X)

    # Check correct init for a given value of m_prior_init
    m_prior_init = rng.rand(n_features)
    bgmm = BayesianGaussianMixture(n_components=n_components,
                                   m_prior_init=m_prior_init).fit(X)
    assert_almost_equal(m_prior_init, bgmm._m_prior)

    # Check correct init for the default value of bem_prior_initta
    bgmm = BayesianGaussianMixture(n_components=n_components).fit(X)
    assert_almost_equal(X.mean(axis=0), bgmm._m_prior)


def test_bayesian_mixture_precisions_prior_initialisation():
    rng = np.random.RandomState(0)
    n_samples, n_features = 10, 2
    X = rng.rand(n_samples, n_features)

    # Check raise message for a bad value of nu_prior_init
    bad_nu_prior_init = n_features - 1.
    bgmm = BayesianGaussianMixture(nu_prior_init=bad_nu_prior_init)
    assert_raise_message(ValueError,
                         "The parameter 'nu_prior_init' should be "
                         "greater than %d, but got %.3f."
                         % (n_features - 1, bad_nu_prior_init),
                         bgmm.fit, X)

    # Check correct init for a given value of nu_prior_init
    nu_prior_init = rng.rand() + n_features - 1.
    bgmm = BayesianGaussianMixture(nu_prior_init=nu_prior_init).fit(X)
    assert_almost_equal(nu_prior_init, bgmm._nu_prior)

    # Check correct init for the default value of nu_prior_init
    nu_prior_init_default = n_features
    bgmm = BayesianGaussianMixture(nu_prior_init=nu_prior_init_default).fit(X)
    assert_almost_equal(nu_prior_init_default, bgmm._nu_prior)

    # Check correct init for a given value of precision_prior_init
    precision_prior_init = {
        'full': np.cov(X.T, bias=1),
        'tied': np.cov(X.T, bias=1),
        'diag': np.diag(np.atleast_2d(np.cov(X.T, bias=1))),
        'spherical': rng.rand()}

    bgmm = BayesianGaussianMixture()
    for cov_type in ['full', 'tied', 'diag', 'spherical']:
        bgmm.covariance_type = cov_type
        bgmm.precision_prior_init = precision_prior_init[cov_type]
        bgmm.fit(X)
        assert_almost_equal(precision_prior_init[cov_type],
                            bgmm._precision_prior)

    # Check raise message for a bad spherical value of precision_prior_init
    bad_precision_init = -1.
    bgmm = BayesianGaussianMixture(covariance_type='spherical',
                                   precision_prior_init=bad_precision_init)
    assert_raise_message(ValueError,
                         "The parameter 'spherical precision_prior_init' "
                         "should be greater than 0., but got %.3f."
                         % bad_precision_init,
                         bgmm.fit, X)

    # Check correct init for the default value of precision_prior_init
    precision_prior_init_default = {
        'full': np.eye(X.shape[1]),
        'tied': np.eye(X.shape[1]),
        'diag': .5 * np.diag(np.atleast_2d(np.cov(X.T, bias=1))),
        'spherical': .5 * np.diag(np.atleast_2d(np.cov(X.T, bias=1))).mean()}

    bgmm = BayesianGaussianMixture()
    for cov_type in ['full', 'tied', 'diag', 'spherical']:
        bgmm.covariance_type = cov_type
        bgmm.fit(X)
        assert_almost_equal(precision_prior_init_default[cov_type],
                            bgmm._precision_prior)


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


def test_bayesian_mixture_means():
    rng = np.random.RandomState(0)
    n_samples, n_features = 10, 2

    X = rng.rand(n_samples, n_features)
    bgmm = BayesianGaussianMixture().fit(X)

    # Check the means values
    assert_almost_equal(bgmm.means_, bgmm.m_)


def test_bayessian_mixture_covariances():
    rng = np.random.RandomState(0)
    n_samples, n_features = 10, 2

    X = rng.rand(n_samples, n_features)
    bgmm = BayesianGaussianMixture().fit(X)

    for covariance_type in ['full', 'tied', 'diag', 'spherical']:
        bgmm.covariance_type = covariance_type
        bgmm.fit(X)

        if covariance_type is 'full':
            pred_covar = bgmm.precisions_ / bgmm.nu_[:, np.newaxis, np.newaxis]
        elif covariance_type is 'diag':
            pred_covar = bgmm.precisions_ / bgmm.nu_[:, np.newaxis]
        else:
            pred_covar = bgmm.precisions_ / bgmm.nu_

        assert_array_almost_equal(pred_covar, bgmm.covariances_)


def generate_data(n_samples, means, covars, random_state=0):
    rng = np.random.RandomState(random_state)
    n_components = len(n_samples)
    X = np.vstack([rng.multivariate_normal(means[j], covars[j], n_samples[j])
                  for j in range(n_components)])
    y = np.concatenate([j * np.ones(n_samples[j])
                       for j in range(n_components)])
    return X, y
