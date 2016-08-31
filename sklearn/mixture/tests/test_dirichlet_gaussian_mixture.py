# Author: Thierry Guillemot <thierry.guillemot.work@gmail.com>
# License: BSD 3 clause

import numpy as np

from sklearn.utils.testing import assert_raise_message
from sklearn.utils.testing import assert_almost_equal

from sklearn.mixture import DirichletGaussianMixture

from sklearn.mixture.tests.test_gaussian_mixture import RandomData
from sklearn.exceptions import ConvergenceWarning
from sklearn.utils.testing import assert_greater_equal, ignore_warnings


COVARIANCE_TYPE = ['full', 'tied', 'diag', 'spherical']


def test_dirichlet_mixture_covariance_type():
    rng = np.random.RandomState(0)
    n_samples, n_features = 10, 2
    X = rng.rand(n_samples, n_features)

    covariance_type = 'bad_covariance_type'
    dpgm = DirichletGaussianMixture(covariance_type=covariance_type,
                                    random_state=rng)
    assert_raise_message(ValueError,
                         "Invalid value for 'covariance_type': %s "
                         "'covariance_type' should be in "
                         "['spherical', 'tied', 'diag', 'full']"
                         % covariance_type,
                         dpgm.fit, X)


def test_dirichlet_mixture_weights_prior_initialisation():
    rng = np.random.RandomState(0)
    n_samples, n_features = 10, 2
    X = rng.rand(n_samples, n_features)

    # Check raise message for a bad value of beta_concentration_prior
    bad_beta_concentration_prior_ = 0.
    dpgm = DirichletGaussianMixture(
        beta_concentration_prior=bad_beta_concentration_prior_,
        random_state=0)
    assert_raise_message(ValueError,
                         "The parameter 'beta_concentration_prior' "
                         "should be greater than 0., but got %.3f."
                         % bad_beta_concentration_prior_,
                         dpgm.fit, X)

    # Check correct init for a given value of beta_concentration_prior
    beta_concentration_prior = rng.rand()
    dpgm = DirichletGaussianMixture(
        beta_concentration_prior=beta_concentration_prior,
        random_state=rng).fit(X)


def test_dirichlet_mixture_means_prior_initialisation():
    rng = np.random.RandomState(0)
    n_samples, n_components, n_features = 10, 3, 2
    X = rng.rand(n_samples, n_features)

    # Check raise message for a bad value of mean_precision_prior
    bad_mean_precision_prior_ = 0.
    dpgm = DirichletGaussianMixture(
        mean_precision_prior=bad_mean_precision_prior_,
        random_state=rng)
    assert_raise_message(ValueError,
                         "The parameter 'mean_precision_prior' should be "
                         "greater than 0., but got %.3f."
                         % bad_mean_precision_prior_,
                         dpgm.fit, X)

    # Check correct init for a given value of mean_precision_prior
    mean_precision_prior = rng.rand()
    dpgm = DirichletGaussianMixture(
        mean_precision_prior=mean_precision_prior,
        random_state=rng).fit(X)
    assert_almost_equal(mean_precision_prior, dpgm.mean_precision_prior_)

    # Check correct init for the default value of mean_precision_prior
    dpgm = DirichletGaussianMixture(random_state=rng).fit(X)
    assert_almost_equal(1., dpgm.mean_precision_prior_)

    # Check raise message for a bad shape of mean_prior
    mean_prior = rng.rand(n_features + 1)
    dpgm = DirichletGaussianMixture(n_components=n_components,
                                    mean_prior=mean_prior,
                                    random_state=rng)
    assert_raise_message(ValueError,
                         "The parameter 'means' should have the shape of ",
                         dpgm.fit, X)

    # Check correct init for a given value of mean_prior
    mean_prior = rng.rand(n_features)
    dpgm = DirichletGaussianMixture(n_components=n_components,
                                    mean_prior=mean_prior,
                                    random_state=rng).fit(X)
    assert_almost_equal(mean_prior, dpgm.mean_prior_)

    # Check correct init for the default value of bemean_priorta
    dpgm = DirichletGaussianMixture(n_components=n_components,
                                    random_state=rng).fit(X)
    assert_almost_equal(X.mean(axis=0), dpgm.mean_prior_)


def test_dirichlet_mixture_precisions_prior_initialisation():
    rng = np.random.RandomState(0)
    n_samples, n_features = 10, 2
    X = rng.rand(n_samples, n_features)

    # Check raise message for a bad value of degrees_of_freedom_prior
    bad_degrees_of_freedom_prior_ = n_features - 1.
    dpgm = DirichletGaussianMixture(
        degrees_of_freedom_prior=bad_degrees_of_freedom_prior_,
        random_state=rng)
    assert_raise_message(ValueError,
                         "The parameter 'degrees_of_freedom_prior' should be "
                         "greater than %d, but got %.3f."
                         % (n_features - 1, bad_degrees_of_freedom_prior_),
                         dpgm.fit, X)

    # Check correct init for a given value of degrees_of_freedom_prior
    degrees_of_freedom_prior = rng.rand() + n_features - 1.
    dpgm = DirichletGaussianMixture(
        degrees_of_freedom_prior=degrees_of_freedom_prior,
        random_state=rng).fit(X)
    assert_almost_equal(degrees_of_freedom_prior,
                        dpgm.degrees_of_freedom_prior_)

    # Check correct init for the default value of degrees_of_freedom_prior
    degrees_of_freedom_prior_default = n_features
    dpgm = DirichletGaussianMixture(
        degrees_of_freedom_prior=degrees_of_freedom_prior_default,
        random_state=rng).fit(X)
    assert_almost_equal(degrees_of_freedom_prior_default,
                        dpgm.degrees_of_freedom_prior_)

    # Check correct init for a given value of covariance_prior
    covariance_prior = {
        'full': np.cov(X.T, bias=1) + 10,
        'tied': np.cov(X.T, bias=1) + 5,
        'diag': np.diag(np.atleast_2d(np.cov(X.T, bias=1))) + 3,
        'spherical': rng.rand()}

    dpgm = DirichletGaussianMixture(random_state=rng)
    for cov_type in ['full', 'tied', 'diag', 'spherical']:
        dpgm.covariance_type = cov_type
        dpgm.covariance_prior = covariance_prior[cov_type]
        dpgm.fit(X)
        assert_almost_equal(covariance_prior[cov_type],
                            dpgm.covariance_prior_)

    # Check raise message for a bad spherical value of covariance_prior
    bad_covariance_prior_ = -1.
    dpgm = DirichletGaussianMixture(covariance_type='spherical',
                                    covariance_prior=bad_covariance_prior_,
                                    random_state=rng)
    assert_raise_message(ValueError,
                         "The parameter 'spherical covariance_prior' "
                         "should be greater than 0., but got %.3f."
                         % bad_covariance_prior_,
                         dpgm.fit, X)

    # Check correct init for the default value of covariance_prior
    covariance_prior_default = {
        'full': np.atleast_2d(np.cov(X.T)),
        'tied': np.atleast_2d(np.cov(X.T)),
        'diag': np.var(X, axis=0, ddof=1),
        'spherical': np.var(X, axis=0, ddof=1).mean()}

    dpgm = DirichletGaussianMixture(random_state=0)
    for cov_type in ['full', 'tied', 'diag', 'spherical']:
        dpgm.covariance_type = cov_type
        dpgm.fit(X)
        assert_almost_equal(covariance_prior_default[cov_type],
                            dpgm.covariance_prior_)


def test_bayesian_mixture_check_is_fitted():
    rng = np.random.RandomState(0)
    n_samples, n_features = 10, 2

    # Check raise message
    dpgm = DirichletGaussianMixture(random_state=rng)
    X = rng.rand(n_samples, n_features)
    assert_raise_message(ValueError,
                         'This DirichletGaussianMixture instance is not '
                         'fitted yet.', dpgm.score, X)


def test_bayesian_mixture_weights():
    rng = np.random.RandomState(0)
    n_samples, n_features = 10, 2

    X = rng.rand(n_samples, n_features)
    dpgm = DirichletGaussianMixture(random_state=rng).fit(X)

    # Check the weights sum = 1
    assert_almost_equal(np.sum(dpgm.weights_), 1.0)


@ignore_warnings(category=ConvergenceWarning)
def test_monotonic_likelihood():
    # We check that each step of the each step of variational inference without
    # regularization improve monotonically the training set of the bound
    rng = np.random.RandomState(0)
    rand_data = RandomData(rng, scale=7)
    n_components = rand_data.n_components

    for covar_type in COVARIANCE_TYPE:
        X = rand_data.X[covar_type]
        dpgm = DirichletGaussianMixture(n_components=2 * n_components,
                                        covariance_type=covar_type,
                                        warm_start=True, max_iter=1,
                                        random_state=rng, tol=1e-4)
        current_lower_bound = -np.infty
        # Do one training iteration at a time so we can make sure that the
        # training log likelihood increases after each iteration.
        for _ in range(600):
            prev_lower_bound = current_lower_bound
            current_lower_bound = dpgm.fit(X).lower_bound_
            assert_greater_equal(current_lower_bound, prev_lower_bound)

            if dpgm.converged_:
                break
        assert(dpgm.converged_)


def test_compare_covar_type():
    # We can compare the 'full' precision with the other cov_type if we apply
    # 1 iter of the M-step (done during _initialize_parameters).
    rng = np.random.RandomState(0)
    rand_data = RandomData(rng, scale=7)
    X = rand_data.X['full']
    n_components = rand_data.n_components
    # Computation of the full_covariance
    dpgm = DirichletGaussianMixture(n_components=2 * n_components,
                                    covariance_type='full',
                                    max_iter=1, random_state=0, tol=1e-7)
    dpgm._check_initial_parameters(X)
    dpgm._initialize_parameters(X)
    full_covariances = (dpgm.covariances_ *
                        dpgm.degrees_of_freedom_[:, np.newaxis, np.newaxis])

    # Check tied_covariance = mean(full_covariances, 0)
    dpgm = DirichletGaussianMixture(n_components=2 * n_components,
                                    covariance_type='tied',
                                    max_iter=1, random_state=0, tol=1e-7)
    dpgm._check_initial_parameters(X)
    dpgm._initialize_parameters(X)

    tied_covariance = dpgm.covariances_ * dpgm.degrees_of_freedom_
    assert_almost_equal(tied_covariance, np.mean(full_covariances, 0))

    # Check diag_covariance = diag(full_covariances)
    dpgm = DirichletGaussianMixture(n_components=2 * n_components,
                                    covariance_type='diag',
                                    max_iter=1, random_state=0, tol=1e-7)
    dpgm._check_initial_parameters(X)
    dpgm._initialize_parameters(X)

    diag_covariances = (dpgm.covariances_ *
                        dpgm.degrees_of_freedom_[:, np.newaxis])
    assert_almost_equal(diag_covariances,
                        np.array([np.diag(cov) for cov in full_covariances]))

    # Check spherical_covariance = np.mean(diag_covariances, 0)
    dpgm = DirichletGaussianMixture(n_components=2 * n_components,
                                    covariance_type='spherical',
                                    max_iter=1, random_state=0, tol=1e-7)
    dpgm._check_initial_parameters(X)
    dpgm._initialize_parameters(X)

    spherical_covariances = dpgm.covariances_ * dpgm.degrees_of_freedom_
    assert_almost_equal(spherical_covariances, np.mean(diag_covariances, 1))


@ignore_warnings(category=ConvergenceWarning)
def test_check_covariance_precision():
    # We check that the dot product of the covariance and the precision
    # matrices is identity.
    rng = np.random.RandomState(0)
    rand_data = RandomData(rng, scale=7)
    n_components, n_features = 2 * rand_data.n_components, 2

    # Computation of the full_covariance
    dpgm = DirichletGaussianMixture(n_components=n_components,
                                    max_iter=100, random_state=rng, tol=1e-3,
                                    reg_covar=0)
    for covar_type in COVARIANCE_TYPE:
        dpgm.covariance_type = covar_type
        dpgm.fit(rand_data.X[covar_type])

        if covar_type == 'full':
            for covar, precision in zip(dpgm.covariances_, dpgm.precisions_):
                assert_almost_equal(np.dot(covar, precision),
                                    np.eye(n_features))
        elif covar_type == 'tied':
            assert_almost_equal(np.dot(dpgm.covariances_, dpgm.precisions_),
                                np.eye(n_features))

        elif covar_type == 'diag':
            assert_almost_equal(dpgm.covariances_ * dpgm.precisions_,
                                np.ones((n_components, n_features)))

        else:
            assert_almost_equal(dpgm.covariances_ * dpgm.precisions_,
                                np.ones(n_components))
