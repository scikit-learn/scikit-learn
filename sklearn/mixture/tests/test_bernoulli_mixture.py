# Author: Wei Xue <xuewei4d@gmail.com>
#         Thierry Guillemot <thierry.guillemot.work@gmail.com>
#         Ben Dilday <ben.dilday@enigma.com>
# License: BSD 3 clause

import warnings
import pytest

import numpy as np

from sklearn.exceptions import ConvergenceWarning
from sklearn.utils.testing import assert_almost_equal
from sklearn.utils.testing import assert_array_equal
from sklearn.utils.testing import assert_equal
from sklearn.utils.testing import assert_greater_equal
from sklearn.utils.testing import assert_raise_message
from sklearn.utils.testing import ignore_warnings

from sklearn.mixture.bernoulli_mixture import (
    _log_no_underflow,
    _estimate_log_bernoulli_prob, _estimate_bernoulli_parameters,
    BernoulliMixture
)


def generate_data(rng, n_samples, weights, means):
    _, n_features = means.shape
    X = []
    for _, (w, m) in enumerate(zip(weights / weights.sum(), means)):
        samples = int(np.round(w * n_samples))
        X.append(rng.binomial(1, m, (samples, n_features)))

    X = np.vstack(X)
    return X


def test_log_no_underflow():
    assert _log_no_underflow(0) != -np.inf
    assert _log_no_underflow(0.123) == np.log(0.123)


def test_estimate_bernoulli_parameters():
    rng = np.random.RandomState(0)
    n_samples, n_features, n_components = 500, 5, 2
    X = rng.binomial(1, 0.5, (n_samples, n_features))
    resp = rng.dirichlet([1] * n_components, n_samples)
    nk, _ = _estimate_bernoulli_parameters(X, resp)

    assert_almost_equal(nk.sum(), n_samples)
    assert_array_equal(nk, resp.mean(0) * n_samples)


def test_estimate_log_bernoulli_prob():
    X = np.array([[0, 1], [1, 0]])
    means = np.array([0.5, 0.5]).reshape((1, 2))
    log_prob = _estimate_log_bernoulli_prob(X, means)
    assert all(log_prob == 2 * np.log(0.5))


def test_bernoulli_mixture_attributes():
    # test bad parameters
    rng = np.random.RandomState(0)
    X = rng.binomial(1, 0.5, (10, 2))

    n_components_bad = 0
    bmm = BernoulliMixture(n_components=n_components_bad)
    assert_raise_message(ValueError,
                         "Invalid value for 'n_components': %d "
                         "Estimation requires at least one component"
                         % n_components_bad, bmm.fit, X)

    tol_bad = -1
    bmm = BernoulliMixture(tol=tol_bad)
    assert_raise_message(ValueError,
                         "Invalid value for 'tol': %.5f "
                         "Tolerance used by the EM must be non-negative"
                         % tol_bad, bmm.fit, X)

    reg_covar_bad = -1
    bmm = BernoulliMixture(reg_covar=reg_covar_bad)
    assert_raise_message(ValueError,
                         "Invalid value for 'reg_covar': %.5f "
                         "regularization on covariance must be "
                         "non-negative" % reg_covar_bad, bmm.fit, X)

    max_iter_bad = 0
    bmm = BernoulliMixture(max_iter=max_iter_bad)
    assert_raise_message(ValueError,
                         "Invalid value for 'max_iter': %d "
                         "Estimation requires at least one iteration"
                         % max_iter_bad, bmm.fit, X)

    n_init_bad = 0
    bmm = BernoulliMixture(n_init=n_init_bad)
    assert_raise_message(ValueError,
                         "Invalid value for 'n_init': %d "
                         "Estimation requires at least one run"
                         % n_init_bad, bmm.fit, X)

    init_params_bad = 'bad_method'
    bmm = BernoulliMixture(init_params=init_params_bad)
    assert_raise_message(ValueError,
                         "Unimplemented initialization method '%s'"
                         % init_params_bad,
                         bmm.fit, X)

    # test good parameters
    n_components, tol, n_init, max_iter, reg_covar = 2, 1e-4, 3, 30, 1e-1
    init_params = 'random'
    bmm = BernoulliMixture(n_components=n_components, tol=tol, n_init=n_init,
                           max_iter=max_iter, reg_covar=reg_covar,
                           init_params=init_params).fit(X)

    assert_equal(bmm.n_components, n_components)
    assert_equal(bmm.tol, tol)
    assert_equal(bmm.reg_covar, reg_covar)
    assert_equal(bmm.max_iter, max_iter)
    assert_equal(bmm.n_init, n_init)
    assert_equal(bmm.init_params, init_params)


def test_check_means():
    rng = np.random.RandomState(0)
    means = np.array([[0.25, 0.5], [0.5, 0.75]])
    weights = np.array([1, 1])
    X = generate_data(rng, 10, weights, means)

    n_components, n_features = means.shape

    g = BernoulliMixture(n_components=n_components)

    # Check means bad shape
    means_bad_shape = rng.rand(n_components + 1, n_features)
    g.means_init = means_bad_shape
    assert_raise_message(ValueError,
                         "The parameter 'means' should have the shape of ",
                         g.fit, X)

    # Check good means matrix
    g.means_init = means
    g.fit(X)
    assert_array_equal(means, g.means_init)


@pytest.mark.filterwarnings('ignore:.*did not converge.*')
@pytest.mark.parametrize("seed", (0, 1, 2))
def test_warm_start(seed):
    random_state = seed
    rng = np.random.RandomState(random_state)
    n_samples, n_features, n_components = 500, 5, 2
    X = rng.binomial(1, 0.5, (n_samples, n_features))

    # Assert the warm_start give the same result for the same number of iter
    g = BernoulliMixture(n_components=n_components, n_init=1, max_iter=2,
                         reg_covar=0, random_state=random_state,
                         warm_start=False, init_params='random')
    h = BernoulliMixture(n_components=n_components, n_init=1, max_iter=1,
                         reg_covar=0, random_state=random_state,
                         warm_start=True, init_params='random')

    g.fit(X)
    score1 = h.fit(X).score(X)
    score2 = h.fit(X).score(X)

    assert_almost_equal(g.weights_, h.weights_)
    assert_almost_equal(g.means_, h.means_)
    assert score2 > score1

    # Assert that by using warm_start we can converge to a good solution
    g = BernoulliMixture(n_components=n_components, n_init=1,
                         max_iter=5, reg_covar=0, random_state=random_state,
                         warm_start=False, tol=1e-6, init_params='random')
    h = BernoulliMixture(n_components=n_components, n_init=1,
                         max_iter=5, reg_covar=0, random_state=random_state,
                         warm_start=True, tol=1e-6, init_params='random')

    g.fit(X)
    assert not g.converged_

    h.fit(X)
    # depending on the data there is large variability in the number of
    # refit necessary to converge due to the complete randomness of the
    # data
    for _ in range(1000):
        h.fit(X)
        if h.converged_:
            break
    assert h.converged_


@ignore_warnings(category=ConvergenceWarning)
def test_convergence_detected_with_warm_start():
    # We check that convergence is detected when warm_start=True
    rng = np.random.RandomState(0)
    n_samples, n_features, n_components = 500, 5, 2
    X = rng.binomial(1, 0.5, (n_samples, n_features))

    for max_iter in (1, 2, 50):
        bmm = BernoulliMixture(n_components=n_components, warm_start=True,
                               max_iter=max_iter, random_state=rng,
                               init_params='random')
        for _ in range(100):
            bmm.fit(X)
            if bmm.converged_:
                break
        assert bmm.converged_
        assert max_iter >= bmm.n_iter_


def test_monotonic_likelihood():
    # We check that each step of the EM without regularization improve
    # monotonically the training set likelihood
    rng = np.random.RandomState(0)
    n_samples, n_features, n_components = 500, 5, 2
    X = rng.binomial(1, 0.5, (n_samples, n_features))

    bmm = BernoulliMixture(n_components=n_components,
                           reg_covar=0,
                           warm_start=True, max_iter=1,
                           random_state=rng,
                           init_params='random',
                           tol=1e-7)
    current_log_likelihood = -np.infty
    with warnings.catch_warnings():
        warnings.simplefilter("ignore", ConvergenceWarning)
        # Do one training iteration at a time so we can make sure that the
        # training log likelihood increases after each iteration.
        for _ in range(1000):
            prev_log_likelihood = current_log_likelihood
            try:
                current_log_likelihood = bmm.fit(X).score(X)
            except ConvergenceWarning:
                pass
            assert_greater_equal(current_log_likelihood,
                                 prev_log_likelihood)

            if bmm.converged_:
                break

        assert bmm.converged_


@ignore_warnings(category=ConvergenceWarning)
def test_init():
    # We check that by increasing the n_init number we have a better solution
    n_samples, n_features, n_components = 500, 5, 2
    for random_state in range(25):
        rng = np.random.RandomState(0)
        X = rng.binomial(1, 0.5, (n_samples, n_features))

        bmm1 = BernoulliMixture(n_components=n_components, n_init=1,
                                max_iter=1, random_state=random_state,
                                init_params='random').fit(X)
        bmm2 = BernoulliMixture(n_components=n_components, n_init=10,
                                max_iter=1, random_state=random_state,
                                init_params='random').fit(X)

        assert bmm2.lower_bound_ >= bmm1.lower_bound_
