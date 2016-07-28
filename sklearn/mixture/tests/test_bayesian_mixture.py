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
from sklearn.mixture.bayesian_mixture import _estimate_wishart_entropy
from sklearn.mixture.bayesian_mixture import _gamma_entropy_spherical
from sklearn.mixture.bayesian_mixture import _gamma_entropy_diag

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

    n_features = 2
    nu = np.abs(rng.rand()) + 1.
    inv_W = linalg.inv(make_spd_matrix(n_features, rng))
    inv_W_chol = linalg.cholesky(inv_W, lower=True)

    expected_norm = (-nu * np.sum(np.log(np.diag(inv_W_chol))) -
                     .5 * n_features * nu * np.log(2.) -
                     .25 * n_features * (n_features - 1) * np.log(np.pi) -
                     np.sum(gammaln(.5 * (nu + 1. -
                                          np.arange(1, n_features + 1.)))))
    predected_norm = _log_wishart_norm(nu, inv_W_chol, n_features)

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
                        _log_wishart_norm(nu, inv_W_chol, n_features))
    predected_entropy = _estimate_wishart_entropy(nu, inv_W_chol, log_lambda,
                                                  n_features)

    assert_almost_equal(expected_entropy, predected_entropy)


def test_gamma_entropy_spherical():
    rng = np.random.RandomState(0)

    n_components = 5
    a = rng.rand(n_components)
    b = rng.rand(n_components)

    expected_entropy = gammaln(a) - (a - 1.) * digamma(a) + np.log(b) + a
    predected_entropy = _gamma_entropy_spherical(a, b)

    assert_almost_equal(expected_entropy, predected_entropy)


def test_gamma_entropy_diag():
    rng = np.random.RandomState(0)

    n_components, n_features = 5, 2
    a = rng.rand(n_components)
    b = rng.rand(n_components, n_features)

    expected_entropy = ((gammaln(a) - (a - 1.) * digamma(a) + a) * len(b) +
                        np.sum(np.log(b)))
    predected_entropy = _gamma_entropy_diag(a, b)

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


# def test_bayesian_mixture_precisions_prior_initialisation():
#     rng = np.random.RandomState(0)
#     n_samples, n_features = 10, 2
#     X = rng.rand(n_samples, n_features)

#     # Check raise message for a bad value of nu_init
#     bad_nu_init = n_features - 1.
#     bgmm = BayesianGaussianMixture(nu_init=bad_nu_init)
#     assert_raise_message(ValueError,
#                          "The parameter 'nu_init' should be "
#                          "greater than %d, but got %.3f."
#                          % (n_features - 1, bad_nu_init),
#                          bgmm.fit, X)

#     # Check correct init for a given value of nu_init
#     nu_init = rng.rand() + n_features - 1.
#     bgmm = BayesianGaussianMixture(nu_init=nu_init).fit(X)
#     assert_almost_equal(nu_init, bgmm._nu_prior)

#     # Check correct init for the default value of nu_init
#     nu_init_default = n_features
#     bgmm = BayesianGaussianMixture(nu_init=nu_init_default).fit(X)
#     assert_almost_equal(nu_init_default, bgmm._nu_prior)

#     # Check correct init for a given value of covariance_init
#     covariance_init = {
#         'full': np.cov(X.T, bias=1),
#         'tied': np.cov(X.T, bias=1),
#         'diag': np.diag(np.atleast_2d(np.cov(X.T, bias=1))),
#         'spherical': rng.rand()}

#     bgmm = BayesianGaussianMixture()
#     for cov_type in ['full', 'tied', 'diag', 'spherical']:
#         print(cov_type)
#         bgmm.covariance_type = cov_type
#         bgmm.covariance_init = covariance_init[cov_type]
#         bgmm.fit(X)
#         assert_almost_equal(covariance_init[cov_type],
#                             bgmm._covariance_prior)

#     # Check raise message for a bad spherical value of covariance_init
#     bad_covariance_init = -1.
#     bgmm = BayesianGaussianMixture(covariance_type='spherical',
#                                    covariance_init=bad_covariance_init)
#     assert_raise_message(ValueError,
#                          "The parameter 'spherical covariance_init' "
#                          "should be greater than 0., but got %.3f."
#                          % bad_covariance_init,
#                          bgmm.fit, X)

#     # Check correct init for the default value of covariance_init
#     covariance_init_default = {
#         'full': np.eye(X.shape[1]),
#         'tied': np.eye(X.shape[1]),
#         'diag': .5 * np.diag(np.atleast_2d(np.cov(X.T, bias=1))),
#         'spherical': .5 * np.diag(np.atleast_2d(np.cov(X.T, bias=1))).mean()}

#     bgmm = BayesianGaussianMixture()
#     for cov_type in ['full', 'tied', 'diag', 'spherical']:
#         bgmm.covariance_type = cov_type
#         bgmm.fit(X)
#         assert_almost_equal(covariance_init_default[cov_type],
#                             bgmm._covariance_prior)


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


def test_check_formula_full():
    # The purpose of that test is to check that the equation are compute
    # correctly by decompositing the fit function step by step.
    # Another function will check that the all process is correctly applied

    # Each blob are well separated
    n_components, n_samples, n_features = 3, 1000, 2
    X, y = make_blobs(centers=n_components, n_samples=n_samples,
                      n_features=n_features, random_state=10)
    resp = np.zeros((n_samples, n_components))
    label = cluster.KMeans(n_clusters=n_components, n_init=1,
                           random_state=0).fit(X).labels_
    resp[np.arange(n_samples), label] = 1

    alpha_init, beta_init, nu_init = 1e0, 5e3, 6e1
    m_init = np.zeros(n_features)
    w_init = np.array([[30, 10.5], [10.5, 9]])
    bgmm = BayesianGaussianMixture(n_components=n_components, max_iter=1,
                                   alpha_init=alpha_init,
                                   beta_init=beta_init,
                                   mean_init=m_init,
                                   covariance_init=w_init,
                                   nu_init=nu_init)
    bgmm._check_initial_parameters(X)
    bgmm._initialize_parameters(X)

    # Check the correctness computation of alpha, beta and nu
    nk, xk, sk = _estimate_gaussian_parameters(
        X, resp, bgmm.reg_covar, 'full')
    alpha_k, beta_k, nu_k = alpha_init + nk, beta_init + nk, nu_init + nk
    assert_almost_equal(alpha_k, bgmm.alpha_)
    assert_almost_equal(beta_k, bgmm.beta_)
    assert_almost_equal(nu_k, bgmm.nu_)

    # Check the correctness means and precision_k
    m_k = ((beta_init * m_init + nk[:, np.newaxis] * xk) /
           beta_k[:, np.newaxis])
    w_k = (w_init + sk * nk[:, np.newaxis, np.newaxis] +
           (beta_init * nk / beta_k)[:, np.newaxis, np.newaxis] *
           np.array([np.outer(x - m_init, x - m_init) for x in xk]))
    assert_almost_equal(m_k, bgmm.means_)
    assert_almost_equal(w_k, bgmm.covariances_)

    # Check the correctness of the resp
    from scipy.special import digamma
    log_pi_k = digamma(alpha_k) - digamma(np.sum(alpha_k))
    log_lambda_k = np.sum((nu_k - np.arange(n_features)))
    assert_almost_equal(log_pi_k, bgmm._estimate_log_weights())



def test_bayesian_mixture_means():
    rng = np.random.RandomState(0)
    n_samples, n_features = 10, 2

    X = rng.rand(n_samples, n_features)
    bgmm = BayesianGaussianMixture().fit(X)

    # Check the means values
    assert_almost_equal(bgmm.means_, bgmm.m_)


# def test_bayessian_mixture_covariances():
#     rng = np.random.RandomState(0)
#     n_samples, n_features = 10, 2

#     X = rng.rand(n_samples, n_features)
#     bgmm = BayesianGaussianMixture().fit(X)

#     for covariance_type in ['full', 'tied', 'diag', 'spherical']:
#         bgmm.covariance_type = covariance_type
#         bgmm.fit(X)

#         if covariance_type is 'full':
#             pred_covar = bgmm.precisions_ / bgmm.nu_[:, np.newaxis, np.newaxis]
#         elif covariance_type is 'diag':
#             pred_covar = bgmm.precisions_ / bgmm.nu_[:, np.newaxis]
#         else:
#             pred_covar = bgmm.precisions_ / bgmm.nu_

#         assert_array_almost_equal(pred_covar, bgmm.covariances_)


def generate_data(n_samples, means, covars, random_state=0):
    rng = np.random.RandomState(random_state)
    n_components = len(n_samples)
    X = np.vstack([rng.multivariate_normal(means[j], covars[j], n_samples[j])
                  for j in range(n_components)])
    y = np.concatenate([j * np.ones(n_samples[j])
                       for j in range(n_components)])
    return X, y
