# Author: Wei Xue <xuewei4d@gmail.com>
#         Thierry Guillemot <thierry.guillemot.work@gmail.com>
# License: BSD 3 clauseimport warnings

import numpy as np
from scipy.special import gammaln

from sklearn.utils.testing import assert_array_almost_equal
from sklearn.utils.testing import assert_raise_message
from sklearn.utils.testing import assert_almost_equal

from sklearn.mixture.bayesian_mixture import _log_dirichlet_norm
from sklearn.mixture.bayesian_mixture import _log_wishart_norm

from sklearn.mixture import BayesianGaussianMixture

from sklearn.mixture.gaussian_mixture import _estimate_gaussian_parameters
from sklearn.mixture.tests.test_gaussian_mixture import RandomData
from sklearn.exceptions import ConvergenceWarning
from sklearn.utils.testing import assert_greater_equal

COVARIANCE_TYPE = ['full', 'tied', 'diag', 'spherical']


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

    # Check raise message for a bad value of alpha_prior
    badalpha_prior_ = 0.
    bgmm = BayesianGaussianMixture(alpha_prior=badalpha_prior_)
    assert_raise_message(ValueError,
                         "The parameter 'alpha_prior' should be "
                         "greater than 0., but got %.3f."
                         % badalpha_prior_,
                         bgmm.fit, X)

    # Check correct init for a given value of alpha_prior
    alpha_prior = rng.rand()
    bgmm = BayesianGaussianMixture(alpha_prior=alpha_prior).fit(X)
    assert_almost_equal(alpha_prior, bgmm.alpha_prior_)

    # Check correct init for the default value of alpha_prior
    bgmm = BayesianGaussianMixture(n_components=n_components).fit(X)
    assert_almost_equal(1. / n_components, bgmm.alpha_prior_)


def test_bayesian_mixture_means_prior_initialisation():
    rng = np.random.RandomState(0)
    n_samples, n_components, n_features = 10, 3, 2
    X = rng.rand(n_samples, n_features)

    # Check raise message for a bad value of beta_prior
    badbeta_prior_ = 0.
    bgmm = BayesianGaussianMixture(beta_prior=badbeta_prior_)
    assert_raise_message(ValueError,
                         "The parameter 'beta_prior' should be "
                         "greater than 0., but got %.3f."
                         % badbeta_prior_,
                         bgmm.fit, X)

    # Check correct init for a given value of beta_prior
    beta_prior = rng.rand()
    bgmm = BayesianGaussianMixture(beta_prior=beta_prior).fit(X)
    assert_almost_equal(beta_prior, bgmm.beta_prior_)

    # Check correct init for the default value of beta_prior
    bgmm = BayesianGaussianMixture().fit(X)
    assert_almost_equal(1., bgmm.beta_prior_)

    # Check raise message for a bad shape of mean_prior
    mean_prior = rng.rand(n_features + 1)
    bgmm = BayesianGaussianMixture(n_components=n_components,
                                   mean_prior=mean_prior)
    assert_raise_message(ValueError,
                         "The parameter 'means' should have the shape of ",
                         bgmm.fit, X)

    # Check correct init for a given value of mean_prior
    mean_prior = rng.rand(n_features)
    bgmm = BayesianGaussianMixture(n_components=n_components,
                                   mean_prior=mean_prior).fit(X)
    assert_almost_equal(mean_prior, bgmm.mean_prior_)

    # Check correct init for the default value of bemean_priorta
    bgmm = BayesianGaussianMixture(n_components=n_components).fit(X)
    assert_almost_equal(X.mean(axis=0), bgmm.mean_prior_)


def test_bayesian_mixture_precisions_prior_initialisation():
    rng = np.random.RandomState(0)
    n_samples, n_features = 10, 2
    X = rng.rand(n_samples, n_features)

    # Check raise message for a bad value of nu_prior
    badnu_prior_ = n_features - 1.
    bgmm = BayesianGaussianMixture(nu_prior=badnu_prior_)
    assert_raise_message(ValueError,
                         "The parameter 'nu_prior' should be "
                         "greater than %d, but got %.3f."
                         % (n_features - 1, badnu_prior_),
                         bgmm.fit, X)

    # Check correct init for a given value of nu_prior
    nu_prior = rng.rand() + n_features - 1.
    bgmm = BayesianGaussianMixture(nu_prior=nu_prior).fit(X)
    assert_almost_equal(nu_prior, bgmm.nu_prior_)

    # Check correct init for the default value of nu_prior
    nu_prior_default = n_features
    bgmm = BayesianGaussianMixture(nu_prior=nu_prior_default).fit(X)
    assert_almost_equal(nu_prior_default, bgmm.nu_prior_)

    # Check correct init for a given value of covariance_prior
    covariance_prior = {
        'full': np.cov(X.T, bias=1) + 10,
        'tied': np.cov(X.T, bias=1) + 5,
        'diag': np.diag(np.atleast_2d(np.cov(X.T, bias=1))) + 3,
        'spherical': rng.rand()}

    bgmm = BayesianGaussianMixture()
    for cov_type in ['full', 'tied', 'diag', 'spherical']:
        print(cov_type)
        bgmm.covariance_type = cov_type
        bgmm.covariance_prior = covariance_prior[cov_type]
        bgmm.fit(X)
        assert_almost_equal(covariance_prior[cov_type],
                            bgmm.covariance_prior_)

    # Check raise message for a bad spherical value of covariance_prior
    badcovariance_prior_ = -1.
    bgmm = BayesianGaussianMixture(covariance_type='spherical',
                                   covariance_prior=badcovariance_prior_)
    assert_raise_message(ValueError,
                         "The parameter 'spherical covariance_prior' "
                         "should be greater than 0., but got %.3f."
                         % badcovariance_prior_,
                         bgmm.fit, X)

    # Check correct init for the default value of covariance_prior
    covariance_prior_default = {
        'full': np.atleast_2d(np.cov(X.T)),
        'tied': np.atleast_2d(np.cov(X.T)),
        'diag': np.var(X, axis=0, ddof=1),
        'spherical': np.var(X, axis=0, ddof=1).mean()}

    bgmm = BayesianGaussianMixture()
    for cov_type in ['full', 'tied', 'diag', 'spherical']:
        bgmm.covariance_type = cov_type
        bgmm.fit(X)
        assert_almost_equal(covariance_prior_default[cov_type],
                            bgmm.covariance_prior_)


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


def test_monotonic_likelihood():
    # We check that each step of the EM without regularization improve
    # monotonically the training set of the bound
    rng = np.random.RandomState(0)
    rand_data = RandomData(rng, scale=7)
    n_components = rand_data.n_components

    for covar_type in COVARIANCE_TYPE:
        X = rand_data.X[covar_type]
        bgmm = BayesianGaussianMixture(n_components=2 * n_components,
                                       covariance_type=covar_type, reg_covar=0,
                                       warm_start=True, max_iter=1,
                                       random_state=rng, tol=1e-7)
        current_lower_bound = -np.infty
        with warnings.catch_warnings():
            warnings.simplefilter("ignore", ConvergenceWarning)
            # Do one training iteration at a time so we can make sure that the
            # training log likelihood increases after each iteration.
            for _ in range(500):
                prev_lower_bound = current_lower_bound
                try:
                    current_lower_bound = bgmm.fit(X).lower_bound_
                except ConvergenceWarning:
                    pass
                assert_greater_equal(current_lower_bound,
                                     prev_lower_bound)

                if bgmm.converged_:
                    break
            assert(bgmm.converged_)


def test_compare_covar_type():
    # We can compare the 'full' precision with the other cov_type if we apply
    # 1 iter of the M-step (done during _initialize_parameters).
    rng = np.random.RandomState(0)
    rand_data = RandomData(rng, scale=7)
    X = rand_data.X['full']
    n_components = rand_data.n_components

    # Computation of the full_covariance
    bgmm = BayesianGaussianMixture(n_components=2 * n_components,
                                   covariance_type='full', reg_covar=0,
                                   max_iter=1, random_state=rng, tol=1e-7)
    bgmm._check_initial_parameters(X)
    bgmm._initialize_parameters(X)
    full_covariances = bgmm.covariances_ * bgmm.nu_[:, np.newaxis, np.newaxis]

    # Check tied_covariance = mean(full_covariances, 0)
    bgmm.covariance_type = 'tied'
    bgmm._check_initial_parameters(X)
    bgmm._initialize_parameters(X)

    tied_covariance = bgmm.covariances_ * bgmm.nu_
    assert_almost_equal(tied_covariance, np.mean(full_covariances, 0))

    # Check diag_covariance = diag(full_covariances)
    bgmm.covariance_type = 'diag'
    bgmm._check_initial_parameters(X)
    bgmm._initialize_parameters(X)

    diag_covariances = bgmm.covariances_ * bgmm.nu_[:, np.newaxis]
    assert_almost_equal(diag_covariances,
                        np.array([np.diag(cov) for cov in full_covariances]))

    # Check spherical_covariance = np.mean(diag_covariances, 0)
    bgmm.covariance_type = 'spherical'
    bgmm._check_initial_parameters(X)
    bgmm._initialize_parameters(X)

    spherical_covariances = bgmm.covariances_ * bgmm.nu_
    assert_almost_equal(spherical_covariances, np.mean(diag_covariances, 1))


def test_check_covariance_precision():
    # We check that the dot product of the covariance and the precision
    # matrices is identity.
    rng = np.random.RandomState(0)
    rand_data = RandomData(rng, scale=7)
    n_components, n_features = 2 * rand_data.n_components, 2

    # Computation of the full_covariance
    bgmm = BayesianGaussianMixture(n_components=n_components, reg_covar=0,
                                   max_iter=100, random_state=rng, tol=1e-3)
    for covar_type in COVARIANCE_TYPE:
        bgmm.covariance_type = covar_type
        bgmm.fit(rand_data.X[covar_type])

        if covar_type == 'full':
            for covar, precision in zip(bgmm.covariances_, bgmm.precisions_):
                assert_almost_equal(np.dot(covar, precision),
                                    np.eye(n_features))
        elif covar_type == 'tied':
            assert_almost_equal(np.dot(bgmm.covariances_, bgmm.precisions_),
                                np.eye(n_features))

        elif covar_type == 'diag':
            assert_almost_equal(bgmm.covariances_ * bgmm.precisions_,
                                np.ones((n_components, n_features)))

        else:
            assert_almost_equal(bgmm.covariances_ * bgmm.precisions_,
                                np.ones(n_components))
