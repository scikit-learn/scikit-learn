import sys

import numpy as np

from scipy import stats

from sklearn.covariance import EmpiricalCovariance
from sklearn.datasets.samples_generator import make_spd_matrix
from sklearn.externals.six.moves import cStringIO as StringIO
from sklearn.metrics.cluster import adjusted_rand_score
from sklearn.mixture.gaussian_mixture import GaussianMixture
from sklearn.mixture.gaussian_mixture import _estimate_gaussian_covariance_diag
from sklearn.mixture.gaussian_mixture import _estimate_gaussian_covariance_full
from sklearn.mixture.gaussian_mixture import (
    _estimate_gaussian_covariance_spherical)
from sklearn.mixture.gaussian_mixture import _estimate_gaussian_covariance_tied
from sklearn.exceptions import ConvergenceWarning, NotFittedError
from sklearn.utils.extmath import fast_logdet
from sklearn.utils.testing import assert_allclose
from sklearn.utils.testing import assert_almost_equal
from sklearn.utils.testing import assert_array_almost_equal
from sklearn.utils.testing import assert_array_equal
from sklearn.utils.testing import assert_equal
from sklearn.utils.testing import assert_greater
from sklearn.utils.testing import assert_raise_message
from sklearn.utils.testing import assert_true
from sklearn.utils.testing import assert_warns_message


COVARIANCE_TYPE = ['full', 'tied', 'diag', 'spherical']


def generate_data(n_samples, n_features, weights, means, covariances,
                  covariance_type):
    rng = np.random.RandomState(0)

    X = []
    if covariance_type == 'spherical':
        for _, (w, m, c) in enumerate(zip(weights, means,
                                          covariances['spherical'])):
            X.append(rng.multivariate_normal(m, c * np.eye(n_features),
                                             int(np.round(w * n_samples))))
    if covariance_type == 'diag':
        for _, (w, m, c) in enumerate(zip(weights, means,
                                          covariances['diag'])):
            X.append(rng.multivariate_normal(m, np.diag(c),
                                             int(np.round(w * n_samples))))
    if covariance_type == 'tied':
        for _, (w, m) in enumerate(zip(weights, means)):
            X.append(rng.multivariate_normal(m, covariances['tied'],
                                             int(np.round(w * n_samples))))
    if covariance_type == 'full':
        for _, (w, m, c) in enumerate(zip(weights, means,
                                          covariances['full'])):
            X.append(rng.multivariate_normal(m, c,
                                             int(np.round(w * n_samples))))

    X = np.vstack(X)
    return X


class RandomData(object):
    def __init__(self, rng, n_samples=500, n_components=2, n_features=2,
                 scale=50):
        self.n_samples = n_samples
        self.n_components = n_components
        self.n_features = n_features

        self.weights = rng.rand(n_components)
        self.weights = self.weights / self.weights.sum()
        self.means = rng.rand(n_components, n_features) * scale
        self.covariances = {
            'spherical': .5 + rng.rand(n_components),
            'diag': (.5 + rng.rand(n_components, n_features)) ** 2,
            'tied': make_spd_matrix(n_features, random_state=rng),
            'full': np.array([make_spd_matrix(
                n_features, random_state=rng) * .5
                for _ in range(n_components)])}

        self.X = dict(zip(COVARIANCE_TYPE, [generate_data(
            n_samples, n_features, self.weights, self.means, self.covariances,
            cov_type) for cov_type in COVARIANCE_TYPE]))
        self.Y = np.hstack([k * np.ones(int(np.round(w * n_samples)))
                            for k, w in enumerate(self.weights)])


def test_gaussian_mixture_parameters():
    # test bad parameters
    rng = np.random.RandomState(0)
    X = rng.rand(10, 2)

    n_init = 0
    gmm = GaussianMixture(n_init=n_init)
    assert_raise_message(ValueError,
                         "Invalid value for 'n_init': %d "
                         "Estimation requires at least one run"
                         % n_init,
                         gmm.fit, X)

    max_iter = 0
    gmm = GaussianMixture(max_iter=max_iter)
    assert_raise_message(ValueError,
                         "Invalid value for 'max_iter': %d "
                         "Estimation requires at least one iteration"
                         % max_iter,
                         gmm.fit, X)

    reg_covar = -1
    gmm = GaussianMixture(reg_covar=reg_covar)
    assert_raise_message(ValueError,
                         "Invalid value for 'reg_covar': %.5f "
                         "regularization on covariance must be "
                         "non-negative" % reg_covar,
                         gmm.fit, X)

    # covariance_type should be in [spherical, diag, tied, full]
    covariance_type = 'bad_covariance_type'
    gmm = GaussianMixture(covariance_type=covariance_type)
    assert_raise_message(ValueError,
                         "Invalid value for 'covariance_type': %s "
                         "'covariance_type' should be in "
                         "['spherical', 'tied', 'diag', 'full']"
                         % covariance_type,
                         gmm.fit, X)

    init_params = 'bad_method'
    gmm = GaussianMixture(init_params=init_params)
    assert_raise_message(ValueError,
                         "Unimplemented initialization method '%s'"
                         % init_params,
                         gmm.fit, X)


def test_check_X():
    from sklearn.mixture.base import _check_X
    rng = np.random.RandomState(0)

    n_samples, n_components, n_features = 10, 2, 2

    X_bad_dim = rng.rand(n_components - 1, n_features)
    assert_raise_message(ValueError,
                         'Expected n_samples >= n_components '
                         'but got n_components = %d, n_samples = %d'
                         % (n_components, X_bad_dim.shape[0]),
                         _check_X, X_bad_dim, n_components)

    X_bad_dim = rng.rand(n_components, n_features + 1)
    assert_raise_message(ValueError,
                         'Expected the input data X have %d features, '
                         'but got %d features'
                         % (n_features, X_bad_dim.shape[1]),
                         _check_X, X_bad_dim, n_components, n_features)

    X = rng.rand(n_samples, n_features)
    assert_array_equal(X, _check_X(X, n_components, n_features))


def test_check_weights():
    rng = np.random.RandomState(0)
    rand_data = RandomData(rng)

    n_components = rand_data.n_components
    X = rand_data.X['full']

    g = GaussianMixture(n_components=n_components)

    # Check bad shape
    weights_bad_shape = rng.rand(n_components, 1)
    g.weights_init = weights_bad_shape
    assert_raise_message(ValueError,
                         "The parameter 'weights' should have the shape of "
                         "(%d,), "
                         "but got %s" % (n_components,
                                         str(weights_bad_shape.shape)),
                         g.fit, X)

    # Check bad range
    weights_bad_range = rng.rand(n_components) + 1
    g.weights_init = weights_bad_range
    assert_raise_message(ValueError,
                         "The parameter 'weights' should be in the range "
                         "[0, 1], but got max value %.5f, min value %.5f"
                         % (np.min(weights_bad_range),
                            np.max(weights_bad_range)),
                         g.fit, X)

    # Check bad normalization
    weights_bad_norm = rng.rand(n_components)
    weights_bad_norm = weights_bad_norm / (weights_bad_norm.sum() + 1)
    g.weights_init = weights_bad_norm
    assert_raise_message(ValueError,
                         "The parameter 'weights' should be normalized, "
                         "but got sum(weights) = %.5f"
                         % np.sum(weights_bad_norm),
                         g.fit, X)

    # Check good weights matrix
    weights = rand_data.weights
    g = GaussianMixture(weights_init=weights, n_components=n_components)
    g.fit(X)
    assert_array_equal(weights, g.weights_init)


def test_check_means():
    rng = np.random.RandomState(0)
    rand_data = RandomData(rng)

    n_components, n_features = rand_data.n_components, rand_data.n_features
    X = rand_data.X['full']

    g = GaussianMixture(n_components=n_components)

    # Check means bad shape
    means_bad_shape = rng.rand(n_components + 1, n_features)
    g.means_init = means_bad_shape
    assert_raise_message(ValueError,
                         "The parameter 'means' should have the shape of ",
                         g.fit, X)

    # Check good means matrix
    means = rand_data.means
    g.means_init = means
    g.fit(X)
    assert_array_equal(means, g.means_init)


def test_check_covariances():
    rng = np.random.RandomState(0)
    rand_data = RandomData(rng)

    n_components, n_features = rand_data.n_components, rand_data.n_features

    # Define the bad covariances for each covariance_type
    covariances_bad_shape = {
        'full': rng.rand(n_components + 1, n_features, n_features),
        'tied': rng.rand(n_features + 1, n_features + 1),
        'diag': rng.rand(n_components + 1, n_features),
        'spherical': rng.rand(n_components + 1)}

    # Define not positive-definite covariances
    covariances_not_pos = rng.rand(n_components, n_features, n_features)
    covariances_not_pos[0] = np.eye(n_features)
    covariances_not_pos[0, 0, 0] = -1.

    covariances_not_positive = {
        'full': covariances_not_pos,
        'tied': covariances_not_pos[0],
        'diag': -1. * np.ones((n_components, n_features)),
        'spherical': -1. * np.ones(n_components)}

    not_positive_errors = {
        'full': 'symmetric, positive-definite',
        'tied': 'symmetric, positive-definite',
        'diag': 'positive',
        'spherical': 'positive'}

    for cov_type in ['full', 'tied', 'diag', 'spherical']:
        X = rand_data.X[cov_type]
        g = GaussianMixture(n_components=n_components,
                            covariance_type=cov_type)

        # Check covariance with bad shapes
        g.covariances_init = covariances_bad_shape[cov_type]
        assert_raise_message(ValueError,
                             "The parameter '%s covariance' should have "
                             "the shape of" % cov_type,
                             g.fit, X)

        # Check not positive covariances
        g.covariances_init = covariances_not_positive[cov_type]
        assert_raise_message(ValueError,
                             "'%s covariance' should be %s"
                             % (cov_type, not_positive_errors[cov_type]),
                             g.fit, X)

        # Check the correct init of covariances_init
        g.covariances_init = rand_data.covariances[cov_type]
        g.fit(X)
        assert_array_equal(rand_data.covariances[cov_type], g.covariances_init)


def test_suffstat_sk_full():
    # compare the EmpiricalCovariance.covariance fitted on X*sqrt(resp)
    # with _sufficient_sk_full, n_components=1
    rng = np.random.RandomState(0)
    n_samples, n_features = 500, 2

    # special case 1, assuming data is "centered"
    X = rng.rand(n_samples, n_features)
    resp = rng.rand(n_samples, 1)
    X_resp = np.sqrt(resp) * X
    nk = np.array([n_samples])
    xk = np.zeros((1, n_features))
    covars_pred = _estimate_gaussian_covariance_full(resp, X, nk, xk, 0)
    ecov = EmpiricalCovariance(assume_centered=True)
    ecov.fit(X_resp)
    assert_almost_equal(ecov.error_norm(covars_pred[0], norm='frobenius'), 0)
    assert_almost_equal(ecov.error_norm(covars_pred[0], norm='spectral'), 0)

    # special case 2, assuming resp are all ones
    resp = np.ones((n_samples, 1))
    nk = np.array([n_samples])
    xk = X.mean().reshape((1, -1))
    covars_pred = _estimate_gaussian_covariance_full(resp, X, nk, xk, 0)
    ecov = EmpiricalCovariance(assume_centered=False)
    ecov.fit(X)
    assert_almost_equal(ecov.error_norm(covars_pred[0], norm='frobenius'), 0)
    assert_almost_equal(ecov.error_norm(covars_pred[0], norm='spectral'), 0)


def test_suffstat_sk_tied():
    # use equation Nk * Sk / N = S_tied
    rng = np.random.RandomState(0)
    n_samples, n_features, n_components = 500, 2, 2

    resp = rng.rand(n_samples, n_components)
    resp = resp / resp.sum(axis=1)[:, np.newaxis]
    X = rng.rand(n_samples, n_features)
    nk = resp.sum(axis=0)
    xk = np.dot(resp.T, X) / nk[:, np.newaxis]
    covars_pred_full = _estimate_gaussian_covariance_full(resp, X, nk, xk, 0)
    covars_pred_full = np.sum(nk[:, np.newaxis, np.newaxis] * covars_pred_full,
                              0) / n_samples

    covars_pred_tied = _estimate_gaussian_covariance_tied(resp, X, nk, xk, 0)
    ecov = EmpiricalCovariance()
    ecov.covariance_ = covars_pred_full
    assert_almost_equal(ecov.error_norm(covars_pred_tied, norm='frobenius'), 0)
    assert_almost_equal(ecov.error_norm(covars_pred_tied, norm='spectral'), 0)


def test_suffstat_sk_diag():
    # test against 'full' case
    rng = np.random.RandomState(0)
    n_samples, n_features, n_components = 500, 2, 2

    resp = rng.rand(n_samples, n_components)
    resp = resp / resp.sum(axis=1)[:, np.newaxis]
    X = rng.rand(n_samples, n_features)
    nk = resp.sum(axis=0)
    xk = np.dot(resp.T, X) / nk[:, np.newaxis]
    covars_pred_full = _estimate_gaussian_covariance_full(resp, X, nk, xk, 0)
    covars_pred_full = np.array([np.diag(np.diag(d)) for d in
                                 covars_pred_full])
    covars_pred_diag = _estimate_gaussian_covariance_diag(resp, X, nk, xk, 0)
    covars_pred_diag = np.array([np.diag(d) for d in covars_pred_diag])
    ecov = EmpiricalCovariance()
    for (cov_full, cov_diag) in zip(covars_pred_full, covars_pred_diag):
        ecov.covariance_ = cov_full
        assert_almost_equal(ecov.error_norm(cov_diag, norm='frobenius'), 0)
        assert_almost_equal(ecov.error_norm(cov_diag, norm='spectral'), 0)


def test_gaussian_suffstat_sk_spherical():
    # computing spherical covariance equals to the variance of one-dimension
    # data after flattening, n_components=1
    rng = np.random.RandomState(0)
    n_samples, n_features = 500, 2

    X = rng.rand(n_samples, n_features)
    X = X - X.mean()
    resp = np.ones((n_samples, 1))
    nk = np.array([n_samples])
    xk = X.mean()
    covars_pred_spherical = _estimate_gaussian_covariance_spherical(resp, X,
                                                                    nk, xk, 0)
    covars_pred_spherical2 = (np.dot(X.flatten().T, X.flatten()) /
                              (n_features * n_samples))
    assert_almost_equal(covars_pred_spherical, covars_pred_spherical2)


def _naive_lmvnpdf_diag(X, means, covars):
    resp = np.empty((len(X), len(means)))
    stds = np.sqrt(covars)
    for i, (mean, std) in enumerate(zip(means, stds)):
        resp[:, i] = stats.norm.logpdf(X, mean, std).sum(axis=1)
    return resp


def test_gaussian_mixture_log_probabilities():
    from sklearn.mixture.gaussian_mixture import (
        _estimate_log_gaussian_prob_full,
        _estimate_log_gaussian_prob_tied,
        _estimate_log_gaussian_prob_diag,
        _estimate_log_gaussian_prob_spherical)

    # test aginst with _naive_lmvnpdf_diag
    rng = np.random.RandomState(0)
    rand_data = RandomData(rng)
    n_samples = 500
    n_features = rand_data.n_features
    n_components = rand_data.n_components

    means = rand_data.means
    covars_diag = rng.rand(n_components, n_features)
    X = rng.rand(n_samples, n_features)
    log_prob_naive = _naive_lmvnpdf_diag(X, means, covars_diag)

    # full covariances
    covars_full = np.array([np.diag(x) for x in covars_diag])

    log_prob = _estimate_log_gaussian_prob_full(X, means, covars_full)
    assert_array_almost_equal(log_prob, log_prob_naive)

    # diag covariances
    log_prob = _estimate_log_gaussian_prob_diag(X, means, covars_diag)
    assert_array_almost_equal(log_prob, log_prob_naive)

    # tied
    covars_tied = covars_full.mean(axis=0)
    log_prob_naive = _naive_lmvnpdf_diag(X, means,
                                         [np.diag(covars_tied)] * n_components)
    log_prob = _estimate_log_gaussian_prob_tied(X, means, covars_tied)
    assert_array_almost_equal(log_prob, log_prob_naive)

    # spherical
    covars_spherical = covars_diag.mean(axis=1)
    log_prob_naive = _naive_lmvnpdf_diag(X, means,
                                         [[k] * n_features for k in
                                          covars_spherical])
    log_prob = _estimate_log_gaussian_prob_spherical(X, means,
                                                     covars_spherical)
    assert_array_almost_equal(log_prob, log_prob_naive)

# skip tests on weighted_log_probabilities, log_weights


def test_gaussian_mixture_estimate_log_prob_resp():
    # test whether responsibilities are normalized
    rng = np.random.RandomState(0)
    rand_data = RandomData(rng)
    n_samples = rand_data.n_samples
    n_features = rand_data.n_features
    n_components = rand_data.n_components

    X = rng.rand(n_samples, n_features)
    for cov_type in COVARIANCE_TYPE:
        weights = rand_data.weights
        means = rand_data.means
        covariances = rand_data.covariances[cov_type]
        g = GaussianMixture(n_components=n_components, random_state=rng,
                            weights_init=weights, means_init=means,
                            covariances_init=covariances,
                            covariance_type=cov_type)
        g._initialize_parameters(X)
        _, _, log_resp = g._estimate_log_prob_resp(X)
        resp = np.exp(log_resp)
        assert_array_almost_equal(resp.sum(axis=1), np.ones(n_samples))


def test_gaussian_mixture_predict_predict_proba():
    rng = np.random.RandomState(0)
    rand_data = RandomData(rng)
    for cov_type in COVARIANCE_TYPE:
        X = rand_data.X[cov_type]
        Y = rand_data.Y
        g = GaussianMixture(n_components=rand_data.n_components,
                            random_state=rng, weights_init=rand_data.weights,
                            means_init=rand_data.means,
                            covariances_init=rand_data.covariances[cov_type],
                            covariance_type=cov_type)

        # Check a warning message arrive if we don't do fit
        assert_raise_message(NotFittedError,
                             "This GaussianMixture instance is not fitted "
                             "yet. Call 'fit' with appropriate arguments "
                             "before using this method.", g.predict, X)

        g.fit(X)
        Y_pred = g.predict(X)
        Y_pred_proba = g.predict_proba(X).argmax(axis=1)
        assert_array_equal(Y_pred, Y_pred_proba)
        assert_greater(adjusted_rand_score(Y, Y_pred), .95)


def test_gaussian_mixture_fit():
    # recover the ground truth
    rng = np.random.RandomState(0)
    rand_data = RandomData(rng)
    n_features = rand_data.n_features
    n_components = rand_data.n_components

    for cov_type in COVARIANCE_TYPE:
        X = rand_data.X[cov_type]
        g = GaussianMixture(n_components=n_components, n_init=20, max_iter=100,
                            reg_covar=0, random_state=rng,
                            covariance_type=cov_type)
        g.fit(X)
        # needs more data to pass the test with rtol=1e-7
        assert_allclose(np.sort(g.weights_), np.sort(rand_data.weights),
                        rtol=0.1, atol=1e-2)

        arg_idx1 = g.means_[:, 0].argsort()
        arg_idx2 = rand_data.means[:, 0].argsort()
        assert_allclose(g.means_[arg_idx1], rand_data.means[arg_idx2],
                        rtol=0.1, atol=1e-2)

        if cov_type == 'spherical':
            cov_pred = np.array([np.eye(n_features) * c
                                 for c in g.covariances_])
            cov_test = np.array([np.eye(n_features) * c for c in
                                 rand_data.covariances['spherical']])
        elif cov_type == 'diag':
            cov_pred = np.array([np.diag(d) for d in g.covariances_])
            cov_test = np.array([np.diag(d) for d in
                                 rand_data.covariances['diag']])
        elif cov_type == 'tied':
            cov_pred = np.array([g.covariances_] * n_components)
            cov_test = np.array([rand_data.covariances['tied']] * n_components)
        elif cov_type == 'full':
            cov_pred = g.covariances_
            cov_test = rand_data.covariances['full']
        arg_idx1 = np.trace(cov_pred, axis1=1, axis2=2).argsort()
        arg_idx2 = np.trace(cov_test, axis1=1, axis2=2).argsort()
        for k, h in zip(arg_idx1, arg_idx2):
            ecov = EmpiricalCovariance()
            ecov.covariance_ = cov_test[h]
            # the accuracy depends on the number of data and randomness, rng
            assert_allclose(ecov.error_norm(cov_pred[k]), 0, atol=0.1)


def test_gaussian_mixture_fit_best_params():
    rng = np.random.RandomState(0)
    rand_data = RandomData(rng)
    n_components = rand_data.n_components
    n_init = 10
    for cov_type in COVARIANCE_TYPE:
        X = rand_data.X[cov_type]
        g = GaussianMixture(n_components=n_components, n_init=1,
                            max_iter=100, reg_covar=0, random_state=rng,
                            covariance_type=cov_type)
        ll = []
        for _ in range(n_init):
            g.fit(X)
            ll.append(g.score(X))
        ll = np.array(ll)
        g_best = GaussianMixture(n_components=n_components,
                                 n_init=n_init, max_iter=100, reg_covar=0,
                                 random_state=rng, covariance_type=cov_type)
        g_best.fit(X)
        assert_almost_equal(ll.min(), g_best.score(X))


def test_gaussian_mixture_fit_convergence_warning():
    rng = np.random.RandomState(0)
    rand_data = RandomData(rng, scale=1)
    n_components = rand_data.n_components
    max_iter = 1
    for cov_type in COVARIANCE_TYPE:
        X = rand_data.X[cov_type]
        g = GaussianMixture(n_components=n_components, n_init=1,
                            max_iter=max_iter, reg_covar=0, random_state=rng,
                            covariance_type=cov_type)
        assert_warns_message(ConvergenceWarning,
                             'Initialization %d did not converged. '
                             'Try different init parameters, '
                             'or increase n_init, tol '
                             'or check for degenerate data.'
                             % max_iter, g.fit, X)


def test_gaussian_mixture_n_parameters():
    # Test that the right number of parameters is estimated
    rng = np.random.RandomState(0)
    n_samples, n_features, n_components = 50, 5, 2
    X = rng.randn(n_samples, n_features)
    n_params = {'spherical': 13, 'diag': 21, 'tied': 26, 'full': 41}
    for cv_type in COVARIANCE_TYPE:
        g = GaussianMixture(
            n_components=n_components, covariance_type=cv_type,
            random_state=rng).fit(X)
        assert_equal(g._n_parameters(), n_params[cv_type])


def test_gaussian_mixture_aic_bic():
    # Test the aic and bic criteria
    rng = np.random.RandomState(0)
    n_samples, n_features, n_components = 50, 3, 2
    X = rng.randn(n_samples, n_features)
    # standard gaussian entropy
    sgh = 0.5 * (fast_logdet(np.cov(X.T, bias=1)) +
                 n_features * (1 + np.log(2 * np.pi)))
    for cv_type in COVARIANCE_TYPE:
        g = GaussianMixture(
            n_components=n_components, covariance_type=cv_type,
            random_state=rng)
        g.fit(X)
        aic = 2 * n_samples * sgh + 2 * g._n_parameters()
        bic = (2 * n_samples * sgh +
               np.log(n_samples) * g._n_parameters())
        bound = n_features / np.sqrt(n_samples)
        assert_true((g.aic(X) - aic) / n_samples < bound)
        assert_true((g.bic(X) - bic) / n_samples < bound)


def test_gaussian_mixture_verbose():
    rng = np.random.RandomState(0)
    rand_data = RandomData(rng)
    n_components = rand_data.n_components
    for cov_type in COVARIANCE_TYPE:
        X = rand_data.X[cov_type]
        g = GaussianMixture(n_components=n_components, n_init=1,
                            max_iter=100, reg_covar=0, random_state=rng,
                            covariance_type=cov_type, verbose=1)
        h = GaussianMixture(n_components=n_components, n_init=1,
                            max_iter=100, reg_covar=0, random_state=rng,
                            covariance_type=cov_type, verbose=2)
        old_stdout = sys.stdout
        sys.stdout = StringIO()
        try:
            g.fit(X)
            h.fit(X)
        finally:
            sys.stdout = old_stdout


def test_warm_start():
    random_state = 0
    rng = np.random.RandomState(random_state)
    n_samples, n_features, n_components = 500, 2, 2
    X = rng.rand(n_samples, n_features)

    # Assert the warm_start give the same result for the same number of iter
    g = GaussianMixture(n_components=n_components, n_init=1,
                        max_iter=2, reg_covar=0, random_state=random_state,
                        warm_start=False)
    g.fit(X)
    h = GaussianMixture(n_components=n_components, n_init=1,
                        max_iter=1, reg_covar=0, random_state=random_state,
                        warm_start=True)
    score1 = h.fit(X).score(X)
    score2 = h.fit(X).score(X)
    print(score1, score2)
    assert_almost_equal(g.weights_, h.weights_)
    assert_almost_equal(g.means_, h.means_)
    assert_almost_equal(g.covariances_, h.covariances_)
    assert_greater(score2, score1)

    # Assert that by using warm_start we can converge to a good solution
    g = GaussianMixture(n_components=n_components, n_init=1,
                        max_iter=5, reg_covar=0, random_state=random_state,
                        warm_start=False)
    g.fit(X)
    h = GaussianMixture(n_components=n_components, n_init=1,
                        max_iter=5, reg_covar=0, random_state=random_state,
                        warm_start=True)
    h.fit(X).fit(X)
    assert_true(not g.converged_)
    assert_true(h.converged_)


def test_score():
    cov_type = 'full'
    rng = np.random.RandomState(0)
    rand_data = RandomData(rng, scale=7)
    n_components = rand_data.n_components
    X = rand_data.X[cov_type]

    # Check the error message if we don't call fit
    gmm1 = GaussianMixture(n_components=n_components, n_init=1,
                           max_iter=1, reg_covar=0, random_state=rng,
                           covariance_type=cov_type)
    assert_raise_message(NotFittedError,
                         "This GaussianMixture instance is not fitted "
                         "yet. Call 'fit' with appropriate arguments "
                         "before using this method.", gmm1.score, X)

    # Check score value
    gmm1.fit(X)
    gmm_score = gmm1.score(X)
    gmm_score_proba = gmm1.score_samples(X).mean()
    assert_almost_equal(gmm_score, gmm_score_proba)

    # Check if the score increase
    gmm2 = GaussianMixture(n_components=n_components, n_init=1,
                           max_iter=1000, reg_covar=0, random_state=rng,
                           covariance_type=cov_type).fit(X)
    assert_greater(gmm2.score(X), gmm1.score(X))


def test_score_samples():
    cov_type = 'full'
    rng = np.random.RandomState(0)
    rand_data = RandomData(rng, scale=7)
    n_components = rand_data.n_components
    X = rand_data.X[cov_type]

    # Check the error message if we don't call fit
    gmm = GaussianMixture(n_components=n_components, n_init=1,
                          max_iter=1, reg_covar=0, random_state=rng,
                          covariance_type=cov_type)
    assert_raise_message(NotFittedError,
                         "This GaussianMixture instance is not fitted "
                         "yet. Call 'fit' with appropriate arguments "
                         "before using this method.", gmm.score_samples, X)

    gmm_score_samples = gmm.fit(X).score_samples(X)
    assert_equal(gmm_score_samples.shape[0], rand_data.n_samples)


def test_monotonic_likelihood():
    rng = np.random.RandomState(0)
    rand_data = RandomData(rng, scale=7)
    n_components = rand_data.n_components

    for cov_type in COVARIANCE_TYPE:
        X = rand_data.X[cov_type]
        gmm = GaussianMixture(n_components=n_components,
                              covariance_type=cov_type, reg_covar=0,
                              warm_start=True, max_iter=20, random_state=rng)
        current_log_likelihood = -np.infty
        for _ in range(100):
            prev_log_likelihood = current_log_likelihood
            current_log_likelihood = gmm.fit(X).score(X)
            assert_greater(current_log_likelihood, prev_log_likelihood)


def test_regularisation():
    rng = np.random.RandomState(0)
    n_samples, n_features = 10, 5

    X = rng.rand(n_samples, n_features)
    X_diag = X.copy()
    X_diag[:, 0] = 0
    X_cov = {'full': X,
             'tied': X,
             'diag': X_diag,
             'spherical': np.zeros((n_samples, n_features))}

    for cov_type in COVARIANCE_TYPE:
        gmm = GaussianMixture(n_components=n_samples, covariance_type=cov_type,
                              reg_covar=0, random_state=rng,)
        assert_raise_message(ValueError,
                             "The algorithm has diverged because of too "
                             "few samples per components. "
                             "Try to decrease the number of components, or "
                             "increase reg_covar.", gmm.fit, X_cov[cov_type])

        gmm.set_params(reg_covar=1e-6).fit(X_cov[cov_type])
