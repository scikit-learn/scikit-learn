import numpy as np
from scipy.special import gamma, digamma

from sklearn.utils.testing import assert_array_almost_equal
from sklearn.utils.testing import assert_array_equal
from sklearn.utils.testing import assert_raise_message
from sklearn.utils.testing import assert_almost_equal
from sklearn.utils.testing import assert_allclose

from sklearn.mixture import DirichletProcessGaussianMixture
from sklearn.mixture.bayesianmixture import _define_prior_shape
from sklearn.mixture.tests.test_gaussianmixture import RandData
from sklearn.covariance import EmpiricalCovariance

rng = np.random.RandomState(0)


def test_log_beta_norm():
    from sklearn.mixture.dpgaussianmixture import _log_beta_norm
    a, b = rng.rand(), rng.rand()
    assert_allclose(np.log(gamma(a + b)/gamma(a)/gamma(b)),
                    _log_beta_norm(a, b))


def test_check_weight_prior():
    param_shape = _define_prior_shape(RandData.n_features, 'full')

    dpgm = DirichletProcessGaussianMixture(n_components=2)
    weight_alpha_bad_shape = rng.rand(1)
    assert_raise_message(
        ValueError, "The parameter 'gamma_prior' should have the shape "
                    "of %s, but got %s"
                    % (param_shape['weight_alpha_prior'],
                       (str(weight_alpha_bad_shape.shape))),
        dpgm._check_weight_prior, weight_alpha_bad_shape,
        param_shape['weight_alpha_prior'])

    weight_alpha_prior_bad_value = -10
    assert_raise_message(
        ValueError,
        "The parameter 'gamma_prior' should be "
        "greater than 0, but got %.3f ."
        % weight_alpha_prior_bad_value,
        dpgm._check_weight_prior, weight_alpha_prior_bad_value,
        param_shape['weight_alpha_prior'])


def test_initialize():
    n_components = RandData.n_components
    dpgm = DirichletProcessGaussianMixture(
        n_components=n_components, precision_type='full', reg_covar=0,
        gamma_prior=0.1, init_params='random', random_state=rng)
    dpgm.n_features = RandData.n_features
    dpgm._initialize(RandData.X['full'])

    dpgm._check_initial_parameters()
    assert dpgm._log_beta_norm_gamma_prior is not None


def test_estimate_weights():
    n_components = 4
    gamma_prior = 0.1
    dpgm = DirichletProcessGaussianMixture(
        n_components=n_components, precision_type='full', reg_covar=0,
        gamma_prior=gamma_prior, init_params='random', random_state=rng)
    nk = rng.rand(n_components) * 100
    gamma_a_pred, gamma_b_pred = dpgm._estimate_weights(nk)
    gamma_a_test = 1 + nk
    gamma_b_test = np.zeros(n_components)
    for k in range(n_components):
        gamma_b_test[k] = gamma_prior
        for i in range(k + 1, n_components):
            gamma_b_test[k] += nk[i]
    assert_array_equal(gamma_a_test, gamma_a_pred)
    assert_array_almost_equal(gamma_b_test, gamma_b_pred)


def test_m_step():
    """Test m-functions on 1-dimension 1-component data"""
    mean = 0
    covar = 1
    n_samples = 100
    X = rng.normal(mean, covar, n_samples)
    X = X.reshape(-1, 1)

    dpgm = DirichletProcessGaussianMixture(
        n_components=1, max_iter=1, reg_covar=0, gamma_prior=0.1,
        precision_type='full', init_params='random', random_state=rng)
    dpgm.n_features = 1
    dpgm._initialize(X)

    assert_array_equal(dpgm.weight_gamma_a_, 1 + n_samples)
    assert_array_equal(dpgm.weight_gamma_b_, 0.1)
    assert dpgm._log_beta_norm_gamma_prior is not None


def test_e_step():
    n_components = 4
    gamma_prior = 0.1
    dpgm = DirichletProcessGaussianMixture(
        n_components=n_components, precision_type='full', reg_covar=0,
        gamma_prior=gamma_prior, init_params='random', random_state=rng)
    nk = rng.rand(n_components) * 100
    dpgm.weight_gamma_a_, dpgm.weight_gamma_b_ = dpgm._estimate_weights(nk)
    gamma_a, gamma_b = dpgm.weight_gamma_a_, dpgm.weight_gamma_b_
    log_weights_pred = dpgm._estimate_log_weights()
    log_weights_test = np.zeros(n_components)
    for k in range(n_components):
        log_weights_test[k] = digamma(gamma_a[k]) - digamma(gamma_a[k] +
                                                            gamma_b[k])
        for i in range(k):
            log_weights_test[k] += (digamma(gamma_b[i]) - digamma(gamma_a[i] +
                                                                  gamma_b[i]))
    assert_almost_equal(log_weights_pred, log_weights_test)


def test_estimate_p_weight():
    mean = 0
    covar = 1
    n_samples = 100
    X = rng.normal(mean, covar, n_samples)
    X = X.reshape(-1, 1)
    gamma_prior = 0.1
    dpgm = DirichletProcessGaussianMixture(
        n_components=1, max_iter=1, reg_covar=0, gamma_prior=gamma_prior,
        precision_type='full', init_params='random', random_state=rng)
    dpgm.n_features = 1
    dpgm._initialize(X)
    _, log_prob, resp, _ = dpgm._estimate_log_prob_resp(X)
    gamma_a, gamma_b = dpgm.weight_gamma_a_, dpgm.weight_gamma_b_

    # _estimate_p_weight
    pweight = dpgm._estimate_p_weight()
    assert_array_equal(pweight,
                       dpgm._log_beta_norm_gamma_prior +
                       (gamma_prior - 1) * (digamma(gamma_b) -
                                            digamma(gamma_a + gamma_b)))

    # _estimate_q_weight
    from sklearn.mixture.dpgaussianmixture import _log_beta_norm
    qweight = dpgm._estimate_q_weight()
    assert_array_equal(qweight,
                       _log_beta_norm(gamma_a, gamma_b) +
                       (gamma_a - 1) * (digamma(gamma_a) -
                                        digamma(gamma_a + gamma_b)) +
                       (gamma_b - 1) * (digamma(gamma_b) -
                                        digamma(gamma_a + gamma_b)))


def test_DPGaussianMixture_fit():
    n_features = RandData.n_features
    n_components = RandData.n_components

    # DPGMM with 'tied' and 'spherical' precision cannot find correct mixtures
    # with 'tied' precision, the 'tied' precision / covariance undermine the
    # ability of selecting the suitable number of components
    # DPGMM with 'spherical' precision cannot pass either.
    for p_type in ['full', 'diag']:
        X = RandData.X[p_type]
        g = DirichletProcessGaussianMixture(
            n_components=n_components * 2, precision_type=p_type, n_init=2,
            max_iter=500, reg_covar=0, gamma_prior=0.1,
            init_params='random', random_state=rng)
        g.fit(X)

        mixture_weight = np.sort(g.weights_)[n_components:]
        assert_allclose(np.sort(RandData.weights), mixture_weight, rtol=0.1)

        argidx_pred = np.argsort(g.weights_)[n_components:]
        argidx_test = np.argsort(RandData.weights)

        assert_allclose(g.means_[argidx_pred],
                        RandData.means[argidx_test], rtol=0.1, atol=1e-2)

        if p_type == 'full':
            cov_pred = g.covars_
            cov_test = RandData.covariances['full']
        elif p_type == 'diag':
            cov_pred = np.array([np.diag(d) for d in g.covars_])
            cov_test = np.array([np.diag(d) for d in
                                 RandData.covariances['diag']])
        for k, h in zip(argidx_pred, argidx_test):
            ecov = EmpiricalCovariance()
            ecov.covariance_ = cov_test[h]
            assert_allclose(ecov.error_norm(cov_pred[k]), 0, atol=0.1)

