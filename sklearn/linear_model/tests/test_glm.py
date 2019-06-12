# Authors: Christian Lorentzen <lorentzen.ch@gmail.com>
#
# License: BSD 3 clause

import numpy as np
from numpy.testing import assert_allclose
import pytest
import scipy as sp
from scipy import linalg, optimize, sparse

from sklearn.datasets import make_classification, make_regression
from sklearn.linear_model import GeneralizedLinearRegressor
from sklearn.linear_model._glm import (
    Link,
    IdentityLink,
    LogLink,
    LogitLink,
    TweedieDistribution,
    NormalDistribution, PoissonDistribution,
    GammaDistribution, InverseGaussianDistribution,
    GeneralizedHyperbolicSecant, BinomialDistribution,
)
from sklearn.linear_model import ElasticNet, LogisticRegression, Ridge
from sklearn.metrics import mean_absolute_error
from sklearn.exceptions import ConvergenceWarning

from sklearn.utils.testing import assert_array_equal


@pytest.fixture(scope="module")
def regression_data():
    X, y = make_regression(n_samples=107,
                           n_features=10,
                           n_informative=80, noise=0.5,
                           random_state=2)
    return X, y


@pytest.mark.parametrize('link', Link.__subclasses__())
def test_link_properties(link):
    """Test link inverse and derivative."""
    rng = np.random.RandomState(42)
    x = rng.rand(100)*100
    link = link()  # instantiate object
    if isinstance(link, LogitLink):
        # careful for large x, note expit(36) = 1
        # limit max eta to 15
        x = x / 100 * 15
    assert_allclose(link.link(link.inverse(x)), x)
    # if f(g(x)) = x, then f'(g(x)) = 1/g'(x)
    assert_allclose(link.derivative(link.inverse(x)),
                    1./link.inverse_derivative(x))

    assert (
      link.inverse_derivative2(x).shape == link.inverse_derivative(x).shape)

    # for LogitLink, in the following x should be between 0 and 1.
    # assert_almost_equal(link.inverse_derivative(link.link(x)),
    #                     1./link.derivative(x), decimal=decimal)


@pytest.mark.parametrize(
    'family, expected',
    [(NormalDistribution(), [True, True, True]),
     (PoissonDistribution(), [False, True, True]),
     (TweedieDistribution(power=1.5), [False, True, True]),
     (GammaDistribution(), [False, False, True]),
     (InverseGaussianDistribution(), [False, False, True]),
     (TweedieDistribution(power=4.5), [False, False, True])])
def test_family_bounds(family, expected):
    """Test the valid range of distributions at -1, 0, 1."""
    result = family.in_y_range([-1, 0, 1])
    assert_array_equal(result, expected)


def test_tweedie_distribution_power():
    with pytest.raises(ValueError, match="no distribution exists"):
        TweedieDistribution(power=0.5)

    with pytest.raises(TypeError, match="must be a real number"):
        TweedieDistribution(power=1j)

    with pytest.raises(TypeError, match="must be a real number"):
        dist = TweedieDistribution()
        dist.power = 1j

    dist = TweedieDistribution()
    assert dist._include_lower_bound is False
    dist.power = 1
    assert dist._include_lower_bound is True


@pytest.mark.parametrize(
    'family, chk_values',
    [(NormalDistribution(), [-1.5, -0.1, 0.1, 2.5]),
     (PoissonDistribution(), [0.1, 1.5]),
     (GammaDistribution(), [0.1, 1.5]),
     (InverseGaussianDistribution(), [0.1, 1.5]),
     (TweedieDistribution(power=-2.5), [0.1, 1.5]),
     (TweedieDistribution(power=-1), [0.1, 1.5]),
     (TweedieDistribution(power=1.5), [0.1, 1.5]),
     (TweedieDistribution(power=2.5), [0.1, 1.5]),
     (TweedieDistribution(power=-4), [0.1, 1.5]),
     (GeneralizedHyperbolicSecant(), [0.1, 1.5])])
def test_deviance_zero(family, chk_values):
    """Test deviance(y,y) = 0 for different families."""
    for x in chk_values:
        assert_allclose(family.deviance(x, x), 0, atol=1e-9)


@pytest.mark.parametrize(
    'family, link',
    [(NormalDistribution(), IdentityLink()),
     (PoissonDistribution(), LogLink()),
     (GammaDistribution(), LogLink()),
     (InverseGaussianDistribution(), LogLink()),
     (TweedieDistribution(power=1.5), LogLink()),
     (TweedieDistribution(power=4.5), LogLink())],
    ids=lambda args: args.__class__.__name__)
def test_fisher_matrix(family, link):
    """Test the Fisher matrix numerically.
    Trick: Use numerical differentiation with y = mu"""
    coef = np.array([-2, 1, 0, 1, 2.5])
    phi = 0.5
    rng = np.random.RandomState(42)
    X = rng.randn(10, 5)
    lin_pred = np.dot(X, coef)
    mu = link.inverse(lin_pred)
    weights = rng.randn(10)**2 + 1
    fisher = family._fisher_matrix(coef=coef, phi=phi, X=X, y=mu,
                                   weights=weights, link=link)
    # check that the Fisher matrix is square and positive definite
    assert fisher.ndim == 2
    assert fisher.shape[0] == fisher.shape[1]
    assert np.all(np.linalg.eigvals(fisher) >= 0)

    approx = np.array([]).reshape(0, coef.shape[0])
    for i in range(coef.shape[0]):
        def f(coef):
            return -family._score(coef=coef, phi=phi, X=X, y=mu,
                                  weights=weights, link=link)[i]
        approx = np.vstack(
            [approx, sp.optimize.approx_fprime(xk=coef, f=f, epsilon=1e-5)])
    assert_allclose(fisher, approx, rtol=1e-3)

    # check the observed information matrix
    oim = family._observed_information(coef=coef, phi=phi, X=X, y=mu,
                                       weights=weights, link=link)
    assert oim.ndim == 2
    assert oim.shape == fisher.shape
    assert_allclose(oim, fisher)


def test_sample_weights_validation():
    """Test the raised errors in the validation of sample_weight."""
    # scalar value but not positive
    X = [[1]]
    y = [1]
    weights = 0
    glm = GeneralizedLinearRegressor(fit_intercept=False)
    with pytest.raises(ValueError):
        glm.fit(X, y, weights)

    # Positive weights are accepted
    glm.fit(X, y, sample_weight=1)

    # 2d array
    weights = [[0]]
    with pytest.raises(ValueError):
        glm.fit(X, y, weights)

    # 1d but wrong length
    weights = [1, 0]
    with pytest.raises(ValueError):
        glm.fit(X, y, weights)

    # 1d but only zeros (sum not greater than 0)
    weights = [0, 0]
    X = [[0], [1]]
    y = [1, 2]
    with pytest.raises(ValueError):
        glm.fit(X, y, weights)

    # 5. 1d but with a negative value
    weights = [2, -1]
    with pytest.raises(ValueError):
        glm.fit(X, y, weights)


@pytest.mark.parametrize('f, fam',
                         [('normal', NormalDistribution()),
                          ('poisson', PoissonDistribution()),
                          ('gamma', GammaDistribution()),
                          ('inverse.gaussian', InverseGaussianDistribution()),
                          ('binomial', BinomialDistribution())])
def test_glm_family_argument(f, fam):
    """Test GLM family argument set as string."""
    y = np.array([0.1, 0.5])  # in range of all distributions
    X = np.array([[1], [2]])
    glm = GeneralizedLinearRegressor(family=f, alpha=0).fit(X, y)
    assert isinstance(glm._family_instance, fam.__class__)

    glm = GeneralizedLinearRegressor(family='not a family',
                                     fit_intercept=False)
    with pytest.raises(ValueError):
        glm.fit(X, y)


@pytest.mark.parametrize('l, link',
                         [('identity', IdentityLink()),
                          ('log', LogLink()),
                          ('logit', LogitLink())])
def test_glm_link_argument(l, link):
    """Test GLM link argument set as string."""
    y = np.array([0.1, 0.5])  # in range of all distributions
    X = np.array([[1], [2]])
    glm = GeneralizedLinearRegressor(family='normal', link=l).fit(X, y)
    assert isinstance(glm._link_instance, link.__class__)

    glm = GeneralizedLinearRegressor(family='normal', link='not a link')
    with pytest.raises(ValueError):
        glm.fit(X, y)


@pytest.mark.parametrize('alpha', ['not a number', -4.2])
def test_glm_alpha_argument(alpha):
    """Test GLM for invalid alpha argument."""
    y = np.array([1, 2])
    X = np.array([[1], [2]])
    glm = GeneralizedLinearRegressor(family='normal', alpha=alpha)
    with pytest.raises(ValueError):
        glm.fit(X, y)


@pytest.mark.parametrize('l1_ratio', ['not a number', -4.2, 1.1, [1]])
def test_glm_l1_ratio_argument(l1_ratio):
    """Test GLM for invalid l1_ratio argument."""
    y = np.array([1, 2])
    X = np.array([[1], [2]])
    glm = GeneralizedLinearRegressor(family='normal', l1_ratio=l1_ratio)
    with pytest.raises(ValueError):
        glm.fit(X, y)


@pytest.mark.parametrize('P1', [['a string', 'a string'], [1, [2]], [1, 2, 3],
                                [-1]])
def test_glm_P1_argument(P1):
    """Test GLM for invalid P1 argument."""
    y = np.array([1, 2])
    X = np.array([[1], [2]])
    glm = GeneralizedLinearRegressor(P1=P1, l1_ratio=0.5, check_input=True)
    with pytest.raises((ValueError, TypeError)):
        glm.fit(X, y)


@pytest.mark.parametrize('P2', ['a string', [1, 2, 3], [[2, 3]],
                                sparse.csr_matrix([1, 2, 3]), [-1]])
def test_glm_P2_argument(P2):
    """Test GLM for invalid P2 argument."""
    y = np.array([1, 2])
    X = np.array([[1], [2]])
    glm = GeneralizedLinearRegressor(P2=P2, check_input=True)
    with pytest.raises(ValueError):
        glm.fit(X, y)


def test_glm_P2_positive_semidefinite():
    """Test GLM for a positive semi-definite P2 argument."""
    n_samples, n_features = 10, 5
    y = np.arange(n_samples)
    X = np.zeros((n_samples, n_features))
    P2 = np.diag([100, 10, 5, 0, -1E-5])
    rng = np.random.RandomState(42)
    # construct random orthogonal matrix Q
    Q, R = linalg.qr(rng.randn(n_features, n_features))
    P2 = Q.T @ P2 @ Q
    glm = GeneralizedLinearRegressor(P2=P2, fit_intercept=False,
                                     check_input=True)
    with pytest.raises(ValueError):
        glm.fit(X, y)

    P2 = sparse.csr_matrix(P2)
    glm = GeneralizedLinearRegressor(P2=P2, fit_intercept=False,
                                     check_input=True)
    with pytest.raises(ValueError):
        glm.fit(X, y)


@pytest.mark.parametrize('fit_intercept', ['not bool', 1, 0, [True]])
def test_glm_fit_intercept_argument(fit_intercept):
    """Test GLM for invalid fit_intercept argument."""
    y = np.array([1, 2])
    X = np.array([[1], [1]])
    glm = GeneralizedLinearRegressor(fit_intercept=fit_intercept)
    with pytest.raises(ValueError):
        glm.fit(X, y)


@pytest.mark.parametrize('solver, l1_ratio',
                         [('not a solver', 0), (1, 0), ([1], 0),
                          ('irls', 0.5), ('lbfgs', 0.5), ('newton-cg', 0.5)])
def test_glm_solver_argument(solver, l1_ratio):
    """Test GLM for invalid solver argument."""
    y = np.array([1, 2])
    X = np.array([[1], [2]])
    glm = GeneralizedLinearRegressor(solver=solver, l1_ratio=l1_ratio)
    with pytest.raises(ValueError):
        glm.fit(X, y)


@pytest.mark.parametrize('max_iter', ['not a number', 0, -1, 5.5, [1]])
def test_glm_max_iter_argument(max_iter):
    """Test GLM for invalid max_iter argument."""
    y = np.array([1, 2])
    X = np.array([[1], [2]])
    glm = GeneralizedLinearRegressor(max_iter=max_iter)
    with pytest.raises(ValueError):
        glm.fit(X, y)


@pytest.mark.parametrize('tol', ['not a number', 0, -1.0, [1e-3]])
def test_glm_tol_argument(tol):
    """Test GLM for invalid tol argument."""
    y = np.array([1, 2])
    X = np.array([[1], [2]])
    glm = GeneralizedLinearRegressor(tol=tol)
    with pytest.raises(ValueError):
        glm.fit(X, y)


@pytest.mark.parametrize('warm_start', ['not bool', 1, 0, [True]])
def test_glm_warm_start_argument(warm_start):
    """Test GLM for invalid warm_start argument."""
    y = np.array([1, 2])
    X = np.array([[1], [1]])
    glm = GeneralizedLinearRegressor(warm_start=warm_start)
    with pytest.raises(ValueError):
        glm.fit(X, y)


@pytest.mark.parametrize('start_params',
                         ['not a start_params', ['zero'], [0, 0, 0],
                          [[0, 0]], ['a', 'b']])
def test_glm_start_params_argument(start_params):
    """Test GLM for invalid start_params argument."""
    y = np.array([1, 2])
    X = np.array([[1], [1]])
    glm = GeneralizedLinearRegressor(start_params=start_params)
    with pytest.raises(ValueError):
        glm.fit(X, y)


@pytest.mark.parametrize('selection', ['not a selection', 1, 0, ['cyclic']])
def test_glm_selection_argument(selection):
    """Test GLM for invalid selection argument"""
    y = np.array([1, 2])
    X = np.array([[1], [1]])
    glm = GeneralizedLinearRegressor(selection=selection)
    with pytest.raises(ValueError):
        glm.fit(X, y)


@pytest.mark.parametrize('random_state', ['a string', 0.5, [0]])
def test_glm_random_state_argument(random_state):
    """Test GLM for invalid random_state argument."""
    y = np.array([1, 2])
    X = np.array([[1], [1]])
    glm = GeneralizedLinearRegressor(random_state=random_state)
    with pytest.raises(ValueError):
        glm.fit(X, y)


@pytest.mark.parametrize('diag_fisher', ['not bool', 1, 0, [True]])
def test_glm_diag_fisher_argument(diag_fisher):
    """Test GLM for invalid diag_fisher arguments."""
    y = np.array([1, 2])
    X = np.array([[1], [1]])
    glm = GeneralizedLinearRegressor(diag_fisher=diag_fisher)
    with pytest.raises(ValueError):
        glm.fit(X, y)


@pytest.mark.parametrize('copy_X', ['not bool', 1, 0, [True]])
def test_glm_copy_X_argument(copy_X):
    """Test GLM for invalid copy_X arguments."""
    y = np.array([1, 2])
    X = np.array([[1], [1]])
    glm = GeneralizedLinearRegressor(copy_X=copy_X)
    with pytest.raises(ValueError):
        glm.fit(X, y)


@pytest.mark.parametrize('check_input', ['not bool', 1, 0, [True]])
def test_glm_check_input_argument(check_input):
    """Test GLM for invalid check_input argument."""
    y = np.array([1, 2])
    X = np.array([[1], [1]])
    glm = GeneralizedLinearRegressor(check_input=check_input)
    with pytest.raises(ValueError):
        glm.fit(X, y)


@pytest.mark.parametrize('solver', ['irls', 'lbfgs', 'newton-cg', 'cd'])
def test_glm_identity_regression(solver):
    """Test GLM regression with identity link on a simple dataset."""
    coef = [1., 2.]
    X = np.array([[1, 1, 1, 1, 1], [0, 1, 2, 3, 4]]).T
    y = np.dot(X, coef)
    glm = GeneralizedLinearRegressor(alpha=0, family='normal', link='identity',
                                     fit_intercept=False, solver=solver,
                                     start_params='zero', tol=1e-7)
    res = glm.fit(X, y)
    assert_allclose(res.coef_, coef, rtol=1e-6)


@pytest.mark.parametrize(
    'family',
    [NormalDistribution(), PoissonDistribution(),
     GammaDistribution(), InverseGaussianDistribution(),
     TweedieDistribution(power=1.5), TweedieDistribution(power=4.5),
     GeneralizedHyperbolicSecant()])
@pytest.mark.parametrize('solver, tol', [('irls', 1e-6),
                                         ('lbfgs', 1e-6),
                                         ('newton-cg', 1e-7),
                                         ('cd', 1e-7)])
def test_glm_log_regression(family, solver, tol):
    """Test GLM regression with log link on a simple dataset."""
    coef = [0.2, -0.1]
    X = np.array([[1, 1, 1, 1, 1], [0, 1, 2, 3, 4]]).T
    y = np.exp(np.dot(X, coef))
    glm = GeneralizedLinearRegressor(
                alpha=0, family=family, link='log', fit_intercept=False,
                solver=solver, start_params='guess', tol=tol)
    res = glm.fit(X, y)
    assert_allclose(res.coef_, coef, rtol=5e-6)


# newton-cg may issue a LineSearchWarning, which we filter out
@pytest.mark.filterwarnings('ignore:The line search algorithm')
@pytest.mark.filterwarnings('ignore:Line Search failed')
@pytest.mark.parametrize('solver, tol', [('irls', 1e-6),
                                         ('lbfgs', 1e-6),
                                         ('newton-cg', 1e-6),
                                         ('cd', 1e-6)])
def test_normal_ridge(solver, tol):
    """Test ridge regression for Normal distributions.

    Compare to test_ridge in test_ridge.py.
    """
    rng = np.random.RandomState(42)
    alpha = 1.0

    # 1. With more samples than features
    n_samples, n_features, n_predict = 100, 7, 10
    X, y, coef = make_regression(n_samples=n_samples+n_predict,
                                 n_features=n_features,
                                 n_informative=n_features-2, noise=0.5,
                                 coef=True, random_state=rng)
    y = y[0:n_samples]
    X, T = X[0:n_samples], X[n_samples:]

    # GLM has 1/(2*n) * Loss + 1/2*L2, Ridge has Loss + L2
    ridge = Ridge(alpha=alpha*n_samples, fit_intercept=True, tol=1e-6,
                  solver='svd', normalize=False)
    ridge.fit(X, y)
    glm = GeneralizedLinearRegressor(alpha=1.0, l1_ratio=0, family='normal',
                                     link='identity', fit_intercept=True,
                                     tol=tol, max_iter=100, solver=solver,
                                     check_input=False, random_state=rng)
    glm.fit(X, y)
    assert glm.coef_.shape == (X.shape[1], )
    assert_allclose(glm.coef_, ridge.coef_, rtol=1e-6)
    assert_allclose(glm.intercept_, ridge.intercept_, rtol=1e-5)
    assert_allclose(glm.predict(T), ridge.predict(T), rtol=1e-6)

    ridge = Ridge(alpha=alpha*n_samples, fit_intercept=False, tol=1e-6,
                  solver='svd', normalize=False)
    ridge.fit(X, y)
    glm = GeneralizedLinearRegressor(alpha=1.0, l1_ratio=0, family='normal',
                                     link='identity', fit_intercept=False,
                                     tol=tol, max_iter=100, solver=solver,
                                     check_input=False, random_state=rng,
                                     fit_dispersion='chisqr')
    glm.fit(X, y)
    assert glm.coef_.shape == (X.shape[1], )
    assert_allclose(glm.coef_, ridge.coef_, rtol=1e-5)
    assert_allclose(glm.intercept_, ridge.intercept_, rtol=1e-6)
    assert_allclose(glm.predict(T), ridge.predict(T), rtol=1e-6)
    mu = glm.predict(X)
    assert_allclose(glm.dispersion_,
                    np.sum((y-mu)**2/(n_samples-n_features)))

    # 2. With more features than samples and sparse
    n_samples, n_features, n_predict = 10, 100, 10
    X, y, coef = make_regression(n_samples=n_samples+n_predict,
                                 n_features=n_features,
                                 n_informative=n_features-2, noise=0.5,
                                 coef=True, random_state=rng)
    y = y[0:n_samples]
    X, T = X[0:n_samples], X[n_samples:]

    # GLM has 1/(2*n) * Loss + 1/2*L2, Ridge has Loss + L2
    ridge = Ridge(alpha=alpha*n_samples, fit_intercept=True, tol=1e-9,
                  solver='sag', normalize=False, max_iter=100000,
                  random_state=42)
    ridge.fit(X, y)
    glm = GeneralizedLinearRegressor(alpha=1.0, l1_ratio=0, family='normal',
                                     link='identity', fit_intercept=True,
                                     tol=tol, max_iter=300, solver=solver,
                                     check_input=False, random_state=rng)
    glm.fit(X, y)
    assert glm.coef_.shape == (X.shape[1], )
    assert_allclose(glm.coef_, ridge.coef_, rtol=5e-6)
    assert_allclose(glm.intercept_, ridge.intercept_, rtol=1e-6)
    assert_allclose(glm.predict(T), ridge.predict(T), rtol=1e-5)

    ridge = Ridge(alpha=alpha*n_samples, fit_intercept=False, tol=1e-7,
                  solver='sag', normalize=False, max_iter=1000,
                  random_state=42)
    ridge.fit(X, y)

    glm = GeneralizedLinearRegressor(alpha=1.0, l1_ratio=0, family='normal',
                                     link='identity', fit_intercept=False,
                                     tol=tol*2, max_iter=300, solver=solver,
                                     check_input=False, random_state=rng)
    glm.fit(X, y)
    assert glm.coef_.shape == (X.shape[1], )
    assert_allclose(glm.coef_, ridge.coef_, rtol=1e-4)
    assert_allclose(glm.intercept_, ridge.intercept_, rtol=1e-5)
    assert_allclose(glm.predict(T), ridge.predict(T), rtol=1e-5)


@pytest.mark.parametrize('solver, tol',
                         [('irls', 1e-7),
                          ('lbfgs', 1e-7),
                          ('newton-cg', 1e-7),
                          ('cd', 1e-7)])
def test_poisson_ridge(solver, tol):
    """Test ridge regression with poisson family and LogLink.

    Compare to R's glmnet"""
    # library("glmnet")
    # options(digits=10)
    # df <- data.frame(a=c(-2,-1,1,2), b=c(0,0,1,1), y=c(0,1,1,2))
    # x <- data.matrix(df[,c("a", "b")])
    # y <- df$y
    # fit <- glmnet(x=x, y=y, alpha=0, intercept=T, family="poisson",
    #               standardize=F, thresh=1e-10, nlambda=10000)
    # coef(fit, s=1)
    # (Intercept) -0.12889386979
    # a            0.29019207995
    # b            0.03741173122
    X = np.array([[-2, -1, 1, 2], [0, 0, 1, 1]]).T
    y = np.array([0, 1, 1, 2])
    rng = np.random.RandomState(42)
    glm = GeneralizedLinearRegressor(alpha=1, l1_ratio=0,
                                     fit_intercept=True, family='poisson',
                                     link='log', tol=tol,
                                     solver=solver, max_iter=300,
                                     random_state=rng)
    glm.fit(X, y)
    assert_allclose(glm.intercept_, -0.12889386979, rtol=1e-5)
    assert_allclose(glm.coef_, [0.29019207995, 0.03741173122], rtol=1e-6)


@pytest.mark.parametrize('diag_fisher', [False, True])
def test_normal_enet(diag_fisher):
    """Test elastic net regression with normal/gaussian family."""
    alpha, l1_ratio = 0.3, 0.7
    n_samples, n_features = 20, 2
    rng = np.random.RandomState(42)
    X = rng.randn(n_samples, n_features).copy(order='F')
    beta = rng.randn(n_features)
    y = 2 + np.dot(X, beta) + rng.randn(n_samples)

    # 1. test normal enet on dense data
    glm = GeneralizedLinearRegressor(alpha=alpha, l1_ratio=l1_ratio,
                                     family='normal', link='identity',
                                     fit_intercept=True, tol=1e-8,
                                     max_iter=100, selection='cyclic',
                                     solver='cd', start_params='zero',
                                     check_input=False,
                                     diag_fisher=diag_fisher)
    glm.fit(X, y)

    enet = ElasticNet(alpha=alpha, l1_ratio=l1_ratio, fit_intercept=True,
                      normalize=False, tol=1e-8, copy_X=True)
    enet.fit(X, y)

    assert_allclose(glm.intercept_, enet.intercept_, rtol=2e-7)
    assert_allclose(glm.coef_, enet.coef_, rtol=5e-5)

    # 2. test normal enet on sparse data
    X = sparse.csc_matrix(X)
    glm.fit(X, y)
    assert_allclose(glm.intercept_, enet.intercept_, rtol=2e-7)
    assert_allclose(glm.coef_, enet.coef_, rtol=5e-5)


def test_poisson_enet():
    """Test elastic net regression with poisson family and LogLink.

    Compare to R's glmnet"""
    # library("glmnet")
    # options(digits=10)
    # df <- data.frame(a=c(-2,-1,1,2), b=c(0,0,1,1), y=c(0,1,1,2))
    # x <- data.matrix(df[,c("a", "b")])
    # y <- df$y
    # fit <- glmnet(x=x, y=y, alpha=0.5, intercept=T, family="poisson",
    #               standardize=F, thresh=1e-10, nlambda=10000)
    # coef(fit, s=1)
    # (Intercept) -0.03550978409
    # a            0.16936423283
    # b            .
    glmnet_intercept = -0.03550978409
    glmnet_coef = [0.16936423283, 0.]
    X = np.array([[-2, -1, 1, 2], [0, 0, 1, 1]]).T
    y = np.array([0, 1, 1, 2])
    rng = np.random.RandomState(42)
    glm = GeneralizedLinearRegressor(alpha=1, l1_ratio=0.5, family='poisson',
                                     link='log', solver='cd', tol=1e-8,
                                     selection='random', random_state=rng,
                                     start_params='guess')
    glm.fit(X, y)
    assert_allclose(glm.intercept_, glmnet_intercept, rtol=2e-6)
    assert_allclose(glm.coef_, glmnet_coef, rtol=2e-7)

    # test results with general optimization procedure
    def obj(coef):
        pd = PoissonDistribution()
        link = LogLink()
        N = y.shape[0]
        mu = link.inverse(X @ coef[1:] + coef[0])
        alpha, l1_ratio = (1, 0.5)
        return 1./(2.*N) * pd.deviance(y, mu) \
            + 0.5 * alpha * (1-l1_ratio) * (coef[1:]**2).sum() \
            + alpha * l1_ratio * np.sum(np.abs(coef[1:]))
    res = optimize.minimize(obj, [0, 0, 0], method='nelder-mead', tol=1e-10,
                            options={'maxiter': 1000, 'disp': False})
    assert_allclose(glm.intercept_, res.x[0], rtol=5e-5)
    assert_allclose(glm.coef_, res.x[1:], rtol=1e-5, atol=1e-9)
    assert_allclose(obj(np.concatenate(([glm.intercept_], glm.coef_))),
                    res.fun, rtol=1e-8)

    # same for start_params='zero' and selection='cyclic'
    # with reduced precision
    glm = GeneralizedLinearRegressor(alpha=1, l1_ratio=0.5, family='poisson',
                                     link='log', solver='cd', tol=1e-5,
                                     selection='cyclic', start_params='zero')
    glm.fit(X, y)
    assert_allclose(glm.intercept_, glmnet_intercept, rtol=1e-4)
    assert_allclose(glm.coef_, glmnet_coef, rtol=1e-4)

    # check warm_start, therefore start with different alpha
    glm = GeneralizedLinearRegressor(alpha=0.005, l1_ratio=0.5,
                                     family='poisson', max_iter=300,
                                     link='log', solver='cd', tol=1e-5,
                                     selection='cyclic', start_params='zero')
    glm.fit(X, y)
    # warm start with original alpha and use of sparse matrices
    glm.warm_start = True
    glm.alpha = 1
    X = sparse.csr_matrix(X)
    glm.fit(X, y)
    assert_allclose(glm.intercept_, glmnet_intercept, rtol=1e-4)
    assert_allclose(glm.coef_, glmnet_coef, rtol=1e-4)


@pytest.mark.parametrize('alpha', [0.01, 0.1, 1, 10])
def test_binomial_enet(alpha):
    """Test elastic net regression with binomial family and LogitLink.

    Compare to LogisticRegression.
    """
    l1_ratio = 0.5
    n_samples = 500
    rng = np.random.RandomState(42)
    X, y = make_classification(n_samples=n_samples, n_classes=2, n_features=6,
                               n_informative=5, n_redundant=0, n_repeated=0,
                               random_state=rng)
    log = LogisticRegression(
        penalty='elasticnet', random_state=rng, fit_intercept=False, tol=1e-6,
        max_iter=1000, l1_ratio=l1_ratio, C=1./(n_samples * alpha),
        solver='saga')
    log.fit(X, y)

    glm = GeneralizedLinearRegressor(
        family=BinomialDistribution(), link=LogitLink(), fit_intercept=False,
        alpha=alpha, l1_ratio=l1_ratio, solver='cd', selection='cyclic',
        tol=1e-7)
    glm.fit(X, y)
    assert_allclose(log.intercept_[0], glm.intercept_, rtol=1e-6)
    assert_allclose(log.coef_[0, :], glm.coef_, rtol=5e-6)


@pytest.mark.parametrize(
        "params",
        [
            {"solver": "irls", "start_params": "guess"},
            {"solver": "irls", "start_params": "zero"},
            {"solver": "lbfgs", "start_params": "guess"},
            {"solver": "lbfgs", "start_params": "zero"},
            {"solver": "newton-cg"},
            {"solver": "cd", "selection": "cyclic", "diag_fisher": False},
            {"solver": "cd", "selection": "cyclic", "diag_fisher": True},
            {"solver": "cd", "selection": "random", "diag_fisher": False},
        ],
        ids=lambda params: ', '.join("%s=%s" % (key, val)
                                     for key,  val in params.items())
)
def test_solver_equivalence(params, regression_data):
    X, y = regression_data
    est_ref = GeneralizedLinearRegressor(random_state=2)
    est_ref.fit(X, y)

    estimator = GeneralizedLinearRegressor(**params)
    estimator.set_params(random_state=2)

    estimator.fit(X, y)

    assert_allclose(estimator.intercept_, est_ref.intercept_, rtol=1e-4)
    assert_allclose(estimator.coef_, est_ref.coef_, rtol=1e-4)
    assert_allclose(
        mean_absolute_error(estimator.predict(X), y),
        mean_absolute_error(est_ref.predict(X), y),
        rtol=1e-4
    )


def test_fit_dispersion(regression_data):
    X, y = regression_data

    est1 = GeneralizedLinearRegressor(random_state=2)
    est1.fit(X, y)
    assert not hasattr(est1, "dispersion_")

    est2 = GeneralizedLinearRegressor(random_state=2, fit_dispersion="chisqr")
    est2.fit(X, y)
    assert isinstance(est2.dispersion_, float)

    est3 = GeneralizedLinearRegressor(
            random_state=2, fit_dispersion="deviance")
    est3.fit(X, y)
    assert isinstance(est3.dispersion_, float)

    assert_allclose(est2.dispersion_,  est3.dispersion_)


@pytest.mark.parametrize("solver", ["irls", "lbfgs", "newton-cg", "cd"])
def test_convergence_warning(solver, regression_data):
    X, y = regression_data

    est = GeneralizedLinearRegressor(solver=solver, random_state=2,
                                     max_iter=1, tol=1e-20)
    with pytest.warns(ConvergenceWarning):
        est.fit(X, y)
