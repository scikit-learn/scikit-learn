import numpy as np
from scipy import linalg
from scipy.special import gamma, digamma

from ...datasets.samples_generator import make_spd_matrix
from ...utils.testing import assert_array_almost_equal
from ...utils.testing import assert_array_equal
from ...utils.testing import assert_raise_message
from ...utils.testing import assert_almost_equal
from ...utils.testing import assert_allclose

from ...mixture import BayesianGaussianMixture
from ...mixture.bayesian_mixture import _define_prior_shape
from ...mixture.tests.test_gaussian_mixture import RandData
from ...covariance import EmpiricalCovariance

rng = np.random.RandomState(0)
PRECISION_TYPE = ['full', 'tied', 'diag', 'spherical']


def test_check_mu_prior():
    from ...mixture.bayesian_mixture import _check_mu_prior

    param_shape = _define_prior_shape(RandData.n_features, 'full')

    # beta should be a scalar, m should be a vector
    mu_beta_prior = rng.rand(*param_shape['mu_beta_prior'])
    mu_m_prior = rng.rand(*param_shape['mu_m_prior'])

    mu_beta_prior_bad_shape = np.array([1])
    assert_raise_message(
        ValueError, "The parameter 'beta_prior' should have the shape of "
                    "%s, but got %s"
                    % (param_shape['mu_beta_prior'],
                       (str(mu_beta_prior_bad_shape.shape))),
        _check_mu_prior, mu_m_prior, mu_beta_prior_bad_shape,
        param_shape['mu_m_prior'], param_shape['mu_beta_prior'])

    mu_beta_prior_bad_value = -1
    assert_raise_message(
        ValueError, "The parameter 'beta_prior' should be greater than 0, "
                    "but got %.3f ."
                    % mu_beta_prior_bad_value,
        _check_mu_prior, mu_m_prior, mu_beta_prior_bad_value,
        param_shape['mu_m_prior'], param_shape['mu_beta_prior'])

    mu_m_prior_bad_shape = rng.rand(RandData.n_features, 2)
    assert_raise_message(
        ValueError, "The parameter 'm_prior' should have the shape of "
                    "%s, but got %s"
                    % (param_shape['mu_m_prior'],
                       (str(mu_m_prior_bad_shape.shape))),
        _check_mu_prior, mu_m_prior_bad_shape, mu_beta_prior,
        param_shape['mu_m_prior'], param_shape['mu_beta_prior'])

    assert_array_equal(
        mu_m_prior,
        _check_mu_prior(mu_m_prior, mu_beta_prior,
                        param_shape['mu_m_prior'],
                        param_shape['mu_beta_prior'])[0]
    )

    assert_array_equal(
        mu_beta_prior,
        _check_mu_prior(mu_m_prior, mu_beta_prior,
                        param_shape['mu_m_prior'],
                        param_shape['mu_beta_prior'])[1]
    )


def test_check_lambda_nu_prior():
    from ...mixture.bayesian_mixture import _check_lambda_nu_prior
    param_shape = _define_prior_shape(RandData.n_features, 'full')

    lambda_nu_prior_bad_shape = np.array([10])
    assert_raise_message(
        ValueError, "The parameter 'nu_prior' should have the shape of "
                    "%s, but got %s"
                    % (param_shape['lambda_nu_prior'],
                       (str(lambda_nu_prior_bad_shape.shape))),
        _check_lambda_nu_prior, lambda_nu_prior_bad_shape,
        param_shape['lambda_nu_prior'], 1)

    lambda_nu_prior_bad_value = 0
    n_features = 2
    assert_raise_message(
        ValueError, "The parameter 'nu_prior' "
                    "should be greater then %d, but got "
                    "value %.3f ."
                    % (n_features - 1, lambda_nu_prior_bad_value),
        _check_lambda_nu_prior, lambda_nu_prior_bad_value,
        param_shape['lambda_nu_prior'], n_features)


def test_check_lambda_W_prior():
    from ...mixture.bayesian_mixture import _check_lambda_prior
    param_shape = _define_prior_shape(RandData.n_features, 'full')

    lambda_W_prior_bad_shape = rng.rand(RandData.n_features + 1,
                                        RandData.n_features + 1)
    lambd_nu_prior = 2
    assert_raise_message(
        ValueError, "The parameter 'inv_W_prior' should have the shape of "
                    "%s, but got %s"
                    % (param_shape['lambda_W_prior'],
                       (str(lambda_W_prior_bad_shape.shape))),
        _check_lambda_prior, lambd_nu_prior, lambda_W_prior_bad_shape,
        param_shape['lambda_nu_prior'], param_shape['lambda_W_prior'], 'full')

    lambda_W_prior_bad_value = np.array([[1, 0], [-1, 0]])
    assert_raise_message(
        ValueError, "The parameter 'inv_W_prior' of Wishart distribution "
                    "should be symmetric, positive-definite.",
        _check_lambda_prior, lambd_nu_prior, lambda_W_prior_bad_value,
        param_shape['lambda_nu_prior'], param_shape['lambda_W_prior'], 'full')

    lambda_W_prior_bad_value = np.array([[1, 1], [1, 0]])
    assert_raise_message(
        ValueError, "The parameter 'inv_W_prior' of Wishart distribution "
                    "should be symmetric, positive-definite.",
        _check_lambda_prior, lambd_nu_prior, lambda_W_prior_bad_value,
        param_shape['lambda_nu_prior'], param_shape['lambda_W_prior'], 'full')

    lambda_W_prior = make_spd_matrix(param_shape['lambda_W_prior'][0], rng)
    assert_array_equal(
        lambda_W_prior,
        _check_lambda_prior(lambd_nu_prior, lambda_W_prior,
                            param_shape['lambda_nu_prior'],
                            param_shape['lambda_W_prior'], 'full')[1])
    assert_array_equal(
        lambd_nu_prior,
        _check_lambda_prior(lambd_nu_prior, lambda_W_prior,
                            param_shape['lambda_nu_prior'],
                            param_shape['lambda_W_prior'], 'full')[0])

    param_shape = _define_prior_shape(RandData.n_features, 'diag')
    lambda_W_prior_bad_shape = rng.rand(RandData.n_features,
                                        RandData.n_features)
    lambda_nu_prior = 10
    assert_raise_message(
        ValueError, "The parameter 'inv_W_prior' should have the shape of "
                    "%s, but got %s"
                    % (param_shape['lambda_W_prior'],
                       lambda_W_prior_bad_shape.shape),
        _check_lambda_prior, lambda_nu_prior, lambda_W_prior_bad_shape,
        param_shape['lambda_nu_prior'], param_shape['lambda_W_prior'], 'diag'
    )

    lambda_W_prior_bad_value = rng.rand(*param_shape['lambda_W_prior'])
    lambda_W_prior_bad_value[0] = -1
    assert_raise_message(
        ValueError, "The parameter 'inv_W_prior' of Gamma distributions "
                    "should be positive.",
        _check_lambda_prior, lambda_nu_prior, lambda_W_prior_bad_value,
        param_shape['lambda_nu_prior'], param_shape['lambda_W_prior'], 'diag'
    )

    lambda_W_prior = np.abs(rng.rand(*param_shape['lambda_W_prior'])) + 1
    assert_array_equal(
        lambda_W_prior,
        _check_lambda_prior(lambda_nu_prior, lambda_W_prior,
                            param_shape['lambda_nu_prior'],
                            param_shape['lambda_W_prior'], 'diag')[1]
    )
    assert_array_equal(
        lambda_nu_prior,
        _check_lambda_prior(lambda_nu_prior, lambda_W_prior,
                            param_shape['lambda_nu_prior'],
                            param_shape['lambda_W_prior'], 'diag')[0]
    )


def test_log_dirichlet_norm():
    from ...mixture.bayesian_mixture import _log_dirichlet_norm
    alpha = rng.rand(2)
    assert_allclose(np.log(gamma(np.sum(alpha)) / np.prod(gamma(alpha))),
                    _log_dirichlet_norm(alpha))


def test_log_wishart_norm():
    from ...mixture.bayesian_mixture import _log_wishart_norm
    n_dim = RandData.n_features
    nu = np.abs(rng.rand()) + 1
    W = make_spd_matrix(n_dim, rng)
    inv_W = linalg.inv(W)
    chol = linalg.cholesky(inv_W, lower=True)
    pred = _log_wishart_norm(n_dim, nu, chol)

    # follow the equation
    sol = (np.power(linalg.det(W), -.5 * nu) /
           (np.power(2, .5 * nu * n_dim) *
            np.power(np.pi, .25 * n_dim * (n_dim - 1)) *
            np.prod([gamma(.5 * (nu + 1 - i)) for i in range(1, n_dim + 1)])))
    # numerical issue causes two numbers are not equal
    assert_almost_equal(np.log(sol), pred)

    log_lambda = rng.rand()


def test_log_gamma_norm_spherical():
    from ...mixture.bayesian_mixture import _log_gamma_norm_spherical
    a = 1 + rng.rand()
    b = 1 + rng.rand()
    assert_array_almost_equal(np.log(np.power(b, -a) / gamma(a)),
                              _log_gamma_norm_spherical(a, 1. / b))


def test_log_gamma_norm_diag():
    from ...mixture.bayesian_mixture import _log_gamma_norm_diag
    a = 1 + rng.rand()
    b = 1 + rng.rand(10)
    assert_array_almost_equal(np.sum(np.log(np.power(1. / b, a) / gamma(a))),
                              _log_gamma_norm_diag(a, 1. / b))


def test_wishart_entropy():
    # it is obvious
    pass


def test_gamma_entropy_spherical():
    # it is obvious
    pass


def test_gamma_entropy_diag():
    # it is obvious
    pass


def test_check_weight_prior():
    param_shape = _define_prior_shape(RandData.n_features, 'full')

    bgm = BayesianGaussianMixture(n_components=2)
    weight_alpha_prior_bad_shape = rng.rand(1)
    assert_raise_message(
        ValueError, "The parameter 'alpha_prior' should have the shape "
                    "of %s, but got %s"
                    % (param_shape['weight_alpha_prior'],
                       (str(weight_alpha_prior_bad_shape.shape))),
        bgm._check_weight_prior, weight_alpha_prior_bad_shape,
        param_shape['weight_alpha_prior'])

    weight_alpha_prior_bad_value = -10
    assert_raise_message(
        ValueError,
        "The parameter 'alpha_prior' should be "
        "greater than 0, but got %.3f ."
        % weight_alpha_prior_bad_value,
        bgm._check_weight_prior, weight_alpha_prior_bad_value,
        param_shape['weight_alpha_prior'])


def test_check_initial_parameters():
    n_features = RandData.n_features
    n_components = RandData.n_components
    param_shape = _define_prior_shape(n_features, 'full')

    alpha_prior = rng.rand()
    m_prior = rng.rand(*param_shape['mu_m_prior'])
    beta_prior = rng.rand(*param_shape['mu_beta_prior'])
    nu_prior = n_features + 1
    inv_W_prior = make_spd_matrix(n_features, rng)

    # cover trivial cases, the inside checking functions are already tested
    bgm = BayesianGaussianMixture(n_components=n_components,
                                  precision_type='full', reg_covar=0,
                                  alpha_prior=alpha_prior, random_state=rng)
    bgm.n_features = n_features
    bgm._check_initial_parameters()

    bgm = BayesianGaussianMixture(n_components=n_components, m_prior=m_prior,
                                  precision_type='full', beta_prior=beta_prior,
                                  random_state=rng)
    bgm.n_features = n_features
    bgm._check_initial_parameters()

    bgm = BayesianGaussianMixture(n_components=n_components, nu_prior=nu_prior,
                                  precision_type='full', reg_covar=0,
                                  inv_W_prior=inv_W_prior, random_state=rng)
    bgm.n_features = n_features
    bgm._check_initial_parameters()

    # test only one of a pair of parameters is set
    bgm = BayesianGaussianMixture(n_components=n_components,
                                  precision_type='full', reg_covar=0,
                                  m_prior=m_prior, random_state=rng)
    assert_raise_message(
        ValueError,
        "The user provided prior parameters 'beta_prior' and 'm_prior' "
        "should both be None or not at the same time.",
        bgm._check_initial_parameters)

    bgm = BayesianGaussianMixture(n_components=n_components,
                                  precision_type='full', reg_covar=0,
                                  nu_prior=nu_prior, random_state=rng)
    assert_raise_message(
        ValueError,
        "The user provided prior parameters 'nu_prior' and 'inv_W_prior' "
        "should both be None or not at the same time.",
        bgm._check_initial_parameters)


def test_estimate_suffstat():
    # it is obvious
    pass


def test_initialize():
    """Test the initial values of parameters are valid"""
    n_components = RandData.n_components

    for p_type in PRECISION_TYPE:
        bgm = BayesianGaussianMixture(n_components=n_components,
                                      precision_type=p_type, reg_covar=0,
                                      init_params='random', random_state=rng)
        bgm.n_features = RandData.n_features
        X = RandData.X[p_type]
        bgm._initialize(X)

        # should not raise any exception for auto initialized parameters
        bgm._check_initial_parameters()

        assert bgm._log_gaussian_norm_prior is not None
        if p_type in ['full', 'tied']:
            assert bgm._log_wishart_norm_prior is not None
        if p_type in ['diag', 'spherical']:
            assert bgm._log_gamma_norm_prior is not None

        # bad X distribution
        X_bad_value = np.array(X)
        X_bad_value[:, 1] = 0
        if p_type in ['full', 'tied']:
            assert_raise_message(
                ValueError,
                "inv_W_prior_ must be symmetric, positive-definite. "
                "Check the data distribution.",
                bgm._initialize, X_bad_value)

        if p_type == 'diag':
            assert_raise_message(
                ValueError,
                "inv_W_prior_ must be greater than 0. "
                "Check the data distribution.",
                bgm._initialize, X_bad_value)

        X_bad_value[:] = 1
        if p_type == 'spherical':
            assert_raise_message(
                ValueError,
                "inv_W_prior_ must be greater than 0. "
                "Check the data distribution.",
                bgm._initialize, X_bad_value)


def test_m_step():
    """Test m-functions on 1-dimension 1-component data"""
    mean = 0
    covar = 1
    n_samples = 100
    X = rng.normal(mean, covar, n_samples)
    X = X.reshape(-1, 1)
    bgm = BayesianGaussianMixture(n_components=1, max_iter=1, reg_covar=0,
                                  precision_type='full', init_params='random',
                                  random_state=rng)

    # set max_iter = 0 to test m_step functions in _initialize_parameters
    bgm.n_features = 1
    bgm._initialize(X)

    # the only one component get all samples
    assert_array_equal(bgm.weight_alpha_, bgm.weight_alpha_prior_ + n_samples)
    assert_array_equal(bgm.mu_beta_, bgm.mu_beta_prior_ + n_samples)
    assert_allclose(bgm.mu_m_,
                    np.atleast_2d((bgm.mu_m_prior_ * bgm.mu_beta_prior_ +
                                   X.sum(0)) / bgm.mu_beta_))
    assert_array_equal(bgm.lambda_nu_, bgm.lambda_nu_prior_ + n_samples)
    assert_array_equal(bgm.lambda_inv_W_.squeeze(), bgm.lambda_inv_W_prior_ +
                       np.cov(X.T, bias=1) * n_samples)

    bgm = BayesianGaussianMixture(n_components=1, max_iter=1, reg_covar=0,
                                  precision_type='tied', init_params='random',
                                  random_state=rng)
    bgm.n_features = 1
    bgm._initialize(X)
    assert_allclose(bgm.lambda_inv_W_, bgm.lambda_inv_W_prior_ +
                    np.cov(X.T, bias=1) * n_samples)

    bgm = BayesianGaussianMixture(n_components=1, max_iter=1, reg_covar=0,
                                  precision_type='diag', init_params='random',
                                  random_state=rng)
    bgm.n_features = 1
    bgm._initialize(X)
    assert_array_equal(bgm.lambda_nu_, bgm.lambda_nu_prior_ + .5 * n_samples)
    assert_allclose(bgm.lambda_inv_W_.squeeze(), bgm.lambda_inv_W_prior_ +
                    .5 * np.cov(X.T, bias=1) * n_samples)

    bgm = BayesianGaussianMixture(n_components=1, max_iter=1, reg_covar=0,
                                  precision_type='spherical',
                                  init_params='random', random_state=rng)
    bgm.n_features = 1
    bgm._initialize(X)
    assert_array_equal(bgm.lambda_nu_, bgm.lambda_nu_prior_ + .5 * n_samples)
    assert_allclose(bgm.lambda_inv_W_.squeeze(), bgm.lambda_inv_W_prior_ +
                    .5 / bgm.n_features * np.cov(X.T, bias=1) * n_samples)


def test_e_step():
    mean = 0
    covar = 1
    n_samples = 100
    X = rng.normal(mean, covar, n_samples)
    X = X.reshape(-1, 1)
    bgm = BayesianGaussianMixture(n_components=1, max_iter=1, reg_covar=0,
                                  precision_type='full', init_params='random',
                                  random_state=rng)
    bgm.n_features = 1
    bgm._initialize(X)
    _, log_prob, resp, _ = bgm._estimate_log_prob_resp(X)
    assert_array_equal(bgm._log_pi, 0.0)
    assert_allclose(bgm._log_lambda, digamma(.5 * bgm.lambda_nu_) +
                    np.log(2) - np.log(bgm.lambda_inv_W_.squeeze()))
    assert_array_equal(bgm._inv_W_chol, np.sqrt(bgm.lambda_inv_W_))
    test_prob = (.5 * bgm._log_lambda - .5 * np.log(2 * np.pi) -
                 .5 * (1 / bgm.mu_beta_ + bgm.lambda_nu_ *
                       (X - bgm.mu_m_) ** 2 / bgm.lambda_inv_W_))
    assert_allclose(log_prob.squeeze(), test_prob.squeeze())
    assert_array_equal(resp.squeeze(), np.ones(n_samples))

    bgm = BayesianGaussianMixture(n_components=1, max_iter=1, reg_covar=0,
                                  precision_type='tied', init_params='random',
                                  random_state=rng)
    bgm.n_features = 1
    bgm._initialize(X)
    _, log_prob, resp, _ = bgm._estimate_log_prob_resp(X)
    assert_allclose(bgm._log_lambda, digamma(.5 * bgm.lambda_nu_) +
                    np.log(2) - np.log(bgm.lambda_inv_W_.squeeze()))
    assert_array_equal(bgm._inv_W_chol, np.sqrt(bgm.lambda_inv_W_))
    # reuse test_prob of 'full'
    assert_allclose(log_prob.squeeze(), test_prob.squeeze())
    assert_array_equal(resp.squeeze(), np.ones(n_samples))

    bgm = BayesianGaussianMixture(n_components=1, max_iter=1, reg_covar=0,
                                  precision_type='diag', init_params='random',
                                  random_state=rng)
    bgm.n_features = 1
    bgm._initialize(X)
    _, log_prob, resp, _ = bgm._estimate_log_prob_resp(X)
    assert_allclose(bgm._log_lambda, digamma(bgm.lambda_nu_) -
                    np.log(bgm.lambda_inv_W_.squeeze()))
    test_prob = (.5 * bgm._log_lambda - .5 * np.log(2 * np.pi) -
                 .5 * (1 / bgm.mu_beta_ + bgm.lambda_nu_ *
                       (X - bgm.mu_m_) ** 2 / bgm.lambda_inv_W_))
    assert_allclose(log_prob.squeeze(), test_prob.squeeze())
    assert_array_equal(resp.squeeze(), np.ones(n_samples))

    bgm = BayesianGaussianMixture(n_components=1, max_iter=1,
                                  precision_type='spherical', reg_covar=0,
                                  init_params='random', random_state=rng)
    bgm.n_features = 1
    bgm._initialize(X)
    _, log_prob, resp, _ = bgm._estimate_log_prob_resp(X)
    assert_allclose(bgm._log_lambda, digamma(bgm.lambda_nu_) -
                    np.log(bgm.lambda_inv_W_.squeeze()))
    # reuse test_prob of 'diag'
    assert_allclose(log_prob.squeeze(), test_prob.squeeze())
    assert_array_equal(resp.squeeze(), np.ones(n_samples))


def test_check_is_fitted():
    # it is obvious
    pass


def test_get_parameters():
    pass


def test_set_parameters():
    pass


def test_estimate_lower_bound():
    mean = 0
    covar = 1
    n_samples = 100
    X = rng.normal(mean, covar, n_samples)
    X = X.reshape(-1, 1)
    bgm = BayesianGaussianMixture(n_components=1, max_iter=1, reg_covar=0,
                                  precision_type='full', init_params='random',
                                  random_state=rng)
    # one m_step and one e_step
    bgm.n_features = 1
    bgm._initialize(X)
    _, log_prob, resp, log_resp = bgm._estimate_log_prob_resp(X)

    # _estimate_p_XZ
    pXZ = bgm._estimate_p_XZ(log_prob, resp)
    test_pXZ = 0.5 * n_samples * (
        bgm._log_lambda - np.log(2 * np.pi) - 1 / bgm.mu_beta_ -
        bgm.lambda_nu_ * np.cov(X.T, bias=1) / bgm.lambda_inv_W_ -
        bgm.lambda_nu_ * (X.mean() - bgm.mu_m_) ** 2 / bgm.lambda_inv_W_)
    # Equation 7.5
    assert_allclose(pXZ, test_pXZ.squeeze())

    # _estimate_p_weight
    pweight = bgm._estimate_p_weight()
    assert_array_equal(pweight, bgm._log_dirichlet_norm_alpha_prior)

    # _estimate_p_mu_lambda
    pml = bgm._estimate_p_mu_lambda()
    test_pml = (.5 * (np.log(bgm.mu_beta_prior_ / (2 * np.pi)) +
                      bgm._log_lambda - bgm.mu_beta_prior_ / bgm.mu_beta_ -
                      bgm.mu_beta_prior_ * bgm.lambda_nu_ *
                      (bgm.mu_m_ - bgm.mu_m_prior_) ** 2 /
                      bgm.lambda_inv_W_) +
                bgm._log_wishart_norm_prior +
                .5 * (bgm.lambda_nu_prior_ - 2) * bgm._log_lambda -
                .5 * bgm.lambda_nu_ * bgm.lambda_inv_W_prior_ /
                bgm.lambda_inv_W_)
    assert_allclose(pml, test_pml)

    # _estimate_q_Z
    qz = bgm._estimate_q_Z(resp, log_resp)
    assert_array_equal(qz, np.zeros((n_samples, 1)))

    # _estimate_q_weight
    qweight = bgm._estimate_q_weight()
    from ...mixture.bayesian_mixture import _log_dirichlet_norm
    assert_array_equal(qweight, _log_dirichlet_norm(bgm.weight_alpha_))

    # _estimate_q_mu_lambda
    qml = bgm._estimate_q_mu_lambda()
    from ...mixture.bayesian_mixture import _wishart_entropy
    test_qml = (.5 * bgm._log_lambda +
                .5 * np.log(bgm.mu_beta_ / (2 * np.pi)) - .5 -
                _wishart_entropy(1, bgm.lambda_nu_[0],
                                 bgm._inv_W_chol[0], bgm._log_lambda[0]))
    assert_array_equal(qml, test_qml)

    bgm = BayesianGaussianMixture(n_components=1, max_iter=1, reg_covar=0,
                                  precision_type='tied', init_params='random',
                                  random_state=rng)
    bgm.n_features = 1
    bgm._initialize(X)
    _, log_prob, resp, log_resp = bgm._estimate_log_prob_resp(X)
    # reuse test_pml of 'full'
    pml = bgm._estimate_p_mu_lambda()
    assert_allclose(pml, test_pml)
    qml = bgm._estimate_q_mu_lambda()
    assert_allclose(qml, test_qml)

    bgm = BayesianGaussianMixture(n_components=1, max_iter=1, reg_covar=0,
                                  precision_type='diag', init_params='random',
                                  random_state=rng)
    bgm.n_features = 1
    bgm._initialize(X)
    _, log_prob, resp, log_resp = bgm._estimate_log_prob_resp(X)
    # reuse test_pml of 'full'
    pml = bgm._estimate_p_mu_lambda()
    test_pml = (.5 * (np.log(bgm.mu_beta_prior_ / (2 * np.pi)) +
                      bgm._log_lambda - bgm.mu_beta_prior_ / bgm.mu_beta_ -
                      bgm.mu_beta_prior_ * bgm.lambda_nu_ *
                      (bgm.mu_m_ - bgm.mu_m_prior_) ** 2 / bgm.lambda_inv_W_) +
                bgm.n_components * bgm._log_gamma_norm_prior +
                (bgm.lambda_nu_prior_ - 1) * bgm._log_lambda -
                bgm.lambda_nu_ * bgm.lambda_inv_W_prior_ / bgm.lambda_inv_W_)
    assert_array_equal(pml, test_pml)
    qml = bgm._estimate_q_mu_lambda()
    from ...mixture.bayesian_mixture import _gamma_entropy_diag
    test_qml = (.5 * bgm._log_lambda +
                .5 * np.log(bgm.mu_beta_ / (2 * np.pi)) -
                .5 -
                _gamma_entropy_diag(bgm.lambda_nu_, bgm.lambda_inv_W_))
    assert_allclose(qml, test_qml)

    bgm = BayesianGaussianMixture(n_components=1, max_iter=1,
                                  precision_type='spherical', reg_covar=0,
                                  init_params='random', random_state=rng)
    bgm.n_features = 1
    bgm._initialize(X)
    _, log_prob, resp, log_resp = bgm._estimate_log_prob_resp(X)
    pml = bgm._estimate_p_mu_lambda()
    # reuse test_pml from 'diag'
    assert_array_equal(pml, test_pml)
    qml = bgm._estimate_q_mu_lambda()
    # reuse test_qml from 'diag'
    assert_allclose(qml, test_qml)


def test_BayesianGaussianMixture_fit():
    n_features = RandData.n_features
    n_components = RandData.n_components

    # BGMM with 'tied' and 'spherical' precision cannot find correct mixtures
    # with 'tied' precision, the 'tied' precision / covariance undermine the
    # ability of selecting the suitable number of components
    # BGMM with 'spherical' precision cannot pass the either
    for p_type in ['full', 'diag']:
        X = RandData.X[p_type]
        g = BayesianGaussianMixture(n_components=n_components * 2,
                                    precision_type=p_type, n_init=2,
                                    max_iter=500, reg_covar=0,
                                    alpha_prior=1e-10, init_params='random',
                                    random_state=rng)
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


def test_GaussianMixture_fit_best_params():
    n_components = RandData.n_components
    n_init = 10
    for p_type in PRECISION_TYPE:
        X = RandData.X[p_type]
        g = BayesianGaussianMixture(n_components=n_components, n_init=1,
                                    max_iter=100, reg_covar=0,
                                    random_state=rng, precision_type=p_type)
        ll = []
        for _ in range(n_init):
            g.fit(X)
            ll.append(g.score(X))
        ll = np.array(ll)
        g_best = BayesianGaussianMixture(n_components=n_components,
                                         n_init=n_init, max_iter=100,
                                         reg_covar=0, random_state=rng,
                                         precision_type=p_type)
        g_best.fit(X)
        assert_almost_equal(ll.min(), g_best.score(X))
