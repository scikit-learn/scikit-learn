import numpy as np
from numpy.testing import assert_allclose
import pytest
import scipy as sp
from scipy import sparse

from sklearn.linear_model.glm import (
    Link,
    IdentityLink,
    LogLink,
    TweedieDistribution,
    NormalDistribution, PoissonDistribution,
    GammaDistribution, InverseGaussianDistribution,
    GeneralizedHyperbolicSecand,
    GeneralizedLinearRegressor)
from sklearn.linear_model import ElasticNet, Ridge

from sklearn.utils.testing import (
    assert_equal, assert_almost_equal,
    assert_array_equal, assert_array_almost_equal,
    assert_raises)


def test_link_properties():
    """Test link inverse and derivative
    """
    rng = np.random.RandomState(0)
    x = rng.rand(100)*100
    # from sklearn.linear_model.glm import Link
    # for link in vars()['Link'].__subclasses__():
    for link in Link.__subclasses__():
        link = link()
        assert_almost_equal(link.link(link.inverse(x)), x, decimal=10)
        assert_almost_equal(link.inverse_derivative(link.link(x)),
                            1/link.derivative(x), decimal=10)


def test_family_bounds():
    """Test the valid range of distributions
    """
    family = NormalDistribution()
    result = family.in_y_range([-1, 0, 1])
    assert_array_equal(result, [True, True, True])

    family = PoissonDistribution()
    result = family.in_y_range([-1, 0, 1])
    assert_array_equal(result, [False, True, True])

    family = TweedieDistribution(power=1.5)
    result = family.in_y_range([-1, 0, 1])
    assert_array_equal(result, [False, True, True])

    family = GammaDistribution()
    result = family.in_y_range([-1, 0, 1])
    assert_array_equal(result, [False, False, True])

    family = InverseGaussianDistribution()
    result = family.in_y_range([-1, 0, 1])
    assert_array_equal(result, [False, False, True])

    family = TweedieDistribution(power=4.5)
    result = family.in_y_range([-1, 0, 1])
    assert_array_equal(result, [False, False, True])


def test_deviance_zero():
    """Test deviance(y,y) = 0 for different families
    """
    for family in [NormalDistribution(), PoissonDistribution(),
                   GammaDistribution(), InverseGaussianDistribution(),
                   TweedieDistribution(power=-2.5),
                   TweedieDistribution(power=-1),
                   TweedieDistribution(power=1.5),
                   TweedieDistribution(power=2.5),
                   TweedieDistribution(power=4),
                   GeneralizedHyperbolicSecand()]:
        assert_almost_equal(family.deviance(0.1, 0.1), 0, decimal=10)
        assert_almost_equal(family.deviance(1.5, 1.5), 0, decimal=10)


def test_fisher_matrix():
    """Test the Fisher matrix numerically.
    Trick: Use numerical differentiation with y = mu"""
    for family in [NormalDistribution(), PoissonDistribution(),
                   GammaDistribution(), InverseGaussianDistribution()]:
        link = LogLink()
        rng = np.random.RandomState(0)
        coef = np.array([-2, 1, 0, 1, 2.5])
        phi = 0.5
        X = rng.randn(10, 5)
        lin_pred = np.dot(X, coef)
        mu = link.inverse(lin_pred)
        weights = rng.randn(10)**2 + 1
        fisher = family._fisher_matrix(coef=coef, phi=phi, X=X, y=mu,
                                       weights=weights, link=link)
        approx = np.array([]).reshape(0, coef.shape[0])
        for i in range(coef.shape[0]):
            def f(coef):
                return -family._score(coef=coef, phi=phi, X=X, y=mu,
                                      weights=weights, link=link)[i]
            approx = np.vstack(
                [approx, sp.optimize.approx_fprime(xk=coef, f=f, epsilon=1e-5)]
                )
        assert_allclose(fisher, approx, rtol=1e-3)


def test_sample_weights_validation():
    """Test the raised errors in the validation of sample_weight"""
    # 1. scalar value but not positive
    X = [[1]]
    y = [1]
    weights = 0
    glm = GeneralizedLinearRegressor(fit_intercept=False)
    assert_raises(ValueError, glm.fit, X, y, weights)

    # 2. 2d array
    weights = [[0]]
    assert_raises(ValueError, glm.fit, X, y, weights)

    # 3. 1d but wrong length
    weights = [1, 0]
    assert_raises(ValueError, glm.fit, X, y, weights)

    # 4. 1d but only zeros (sum not greater than 0)
    weights = [0, 0]
    X = [[0], [1]]
    y = [1, 2]
    assert_raises(ValueError, glm.fit, X, y, weights)

    # 5. 1d but weith a negative value
    weights = [2, -1]
    assert_raises(ValueError, glm.fit, X, y, weights)


def test_glm_family_argument():
    """Test GLM family argument set as string
    """
    y = np.array([1, 2])
    X = np.array([[1], [1]])
    for (f, fam) in [('normal', NormalDistribution()),
                     ('poisson', PoissonDistribution()),
                     ('gamma', GammaDistribution()),
                     ('inverse.gaussian', InverseGaussianDistribution())]:
        glm = GeneralizedLinearRegressor(family=f, fit_intercept=False,
                                         alpha=0).fit(X, y)
        assert_equal(type(glm._family_instance), type(fam))

    glm = GeneralizedLinearRegressor(family='not a family',
                                     fit_intercept=False)
    assert_raises(ValueError, glm.fit, X, y)


def test_glm_link_argument():
    """Test GLM link argument set as string
    """
    y = np.array([1, 2])
    X = np.array([[1], [1]])
    for (l, link) in [('identity', IdentityLink()),
                      ('log', LogLink())]:
        glm = GeneralizedLinearRegressor(family='normal', fit_intercept=False,
                                         link=l).fit(X, y)
        assert_equal(type(glm._link_instance), type(link))

    glm = GeneralizedLinearRegressor(family='normal', fit_intercept=False,
                                     link='not a link')
    assert_raises(ValueError, glm.fit, X, y)


def test_glm_alpha_argument():
    """Test GLM alpha argument
    """
    y = np.array([1, 2])
    X = np.array([[1], [1]])
    for alpha in ['not a number', -4.2]:
        glm = GeneralizedLinearRegressor(family='normal', fit_intercept=False,
                                         alpha=alpha)
        assert_raises(ValueError, glm.fit, X, y)


def test_glm_l1_ratio_argument():
    """Test GLM l1_ratio argument
    """
    y = np.array([1, 2])
    X = np.array([[1], [1]])
    for l1_ratio in ['not a number', -4.2, 1.1, [1]]:
        glm = GeneralizedLinearRegressor(family='normal', fit_intercept=False,
                                         l1_ratio=l1_ratio)
        assert_raises(ValueError, glm.fit, X, y)


@pytest.mark.parametrize('P1', [['a string', 'a string'], [1, [2]], [1, 2, 3]])
def test_glm_P1_argument(P1):
    """Test GLM P1 arguments."""
    y = np.array([1, 2])
    X = np.array([[1], [1]])
    glm = GeneralizedLinearRegressor(P1=P1)
    with pytest.raises((ValueError, TypeError)):
        glm.fit(X, y)


@pytest.mark.parametrize('P2', ['a string', [1, 2, 3], [[2, 3]],
                                sparse.csr_matrix([1, 2, 3]),
                                sparse.lil_matrix([[1]])])
def test_glm_P2_argument(P2):
    """Test GLM P2 arguments."""
    y = np.array([1, 2])
    X = np.array([[1], [1]])
    glm = GeneralizedLinearRegressor(P2=P2, fit_intercept=False)
    with pytest.raises((ValueError, TypeError)):
        glm.fit(X, y)


def test_glm_fit_intercept_argument():
    """Test GLM fit_intercept argument
    """
    y = np.array([1, 2])
    X = np.array([[1], [1]])
    for fit_intercept in ['not bool', 1, 0, [True]]:
        glm = GeneralizedLinearRegressor(fit_intercept=fit_intercept)
        assert_raises(ValueError, glm.fit, X, y)


def test_glm_solver_argument():
    """Test GLM solver argument
    """
    y = np.array([1, 2])
    X = np.array([[1], [1]])
    for solver in ['not a solver', 1, [1]]:
        glm = GeneralizedLinearRegressor(solver=solver)
        assert_raises(ValueError, glm.fit, X, y)

    # solver not suitable for L1 penalty
    for solver in ['irls', 'lbfgs', 'newton-cg']:
        glm = GeneralizedLinearRegressor(solver=solver, alpha=1, l1_ratio=0.1)
        assert_raises(ValueError, glm.fit, X, y)


def test_glm_max_iter_argument():
    """Test GLM max_iter argument
    """
    y = np.array([1, 2])
    X = np.array([[1], [1]])
    for max_iter in ['not a number', 0, -1, 5.5, [1]]:
        glm = GeneralizedLinearRegressor(max_iter=max_iter)
        assert_raises(ValueError, glm.fit, X, y)


def test_glm_tol_argument():
    """Test GLM tol argument
    """
    y = np.array([1, 2])
    X = np.array([[1], [1]])
    for tol in ['not a number', 0, -1.0, [1e-3]]:
        glm = GeneralizedLinearRegressor(tol=tol)
        assert_raises(ValueError, glm.fit, X, y)


def test_glm_warm_start_argument():
    """Test GLM warm_start argument
    """
    y = np.array([1, 2])
    X = np.array([[1], [1]])
    for warm_start in ['not bool', 1, 0, [True]]:
        glm = GeneralizedLinearRegressor(warm_start=warm_start)
        assert_raises(ValueError, glm.fit, X, y)


def test_glm_start_params_argument():
    """Test GLM start_params argument
    """
    y = np.array([1, 2])
    X = np.array([[1], [1]])
    for start_params in ['not a start_params', ['zero'], [0, 0, 0],
                         [[0, 0]], ['a', 'b']]:
        glm = GeneralizedLinearRegressor(start_params=start_params)
        assert_raises(ValueError, glm.fit, X, y)


def test_glm_selection_argument():
    """Test GLM selection argument
    """
    y = np.array([1, 2])
    X = np.array([[1], [1]])
    for selection in ['not a selection', 1, 0, ['cyclic']]:
        glm = GeneralizedLinearRegressor(selection=selection)
        assert_raises(ValueError, glm.fit, X, y)


def test_glm_random_state_argument():
    """Test GLM random_state argument
    """
    y = np.array([1, 2])
    X = np.array([[1], [1]])
    for random_state in ['a string', 0.5, [0]]:
        glm = GeneralizedLinearRegressor(random_state=random_state)
        assert_raises(ValueError, glm.fit, X, y)


def test_glm_copy_X_argument():
    """Test GLM copy_X arguments
    """
    y = np.array([1, 2])
    X = np.array([[1], [1]])
    for copy_X in ['not bool', 1, 0, [True]]:
        glm = GeneralizedLinearRegressor(copy_X=copy_X)
        assert_raises(ValueError, glm.fit, X, y)


def test_glm_check_input_argument():
    """Test GLM check_input argument
    """
    y = np.array([1, 2])
    X = np.array([[1], [1]])
    for check_input in ['not bool', 1, 0, [True]]:
        glm = GeneralizedLinearRegressor(check_input=check_input)
        assert_raises(ValueError, glm.fit, X, y)


# TODO: check additional validations if check_input == True

def test_glm_identiy_regression():
    """Test GLM regression with identity link on a simple dataset
    """
    coef = [1, 2]
    X = np.array([[1, 1, 1, 1, 1], [0, 1, 2, 3, 4]]).T
    y = np.dot(X, coef)
    families = (
        NormalDistribution(), PoissonDistribution(),
        GammaDistribution(), InverseGaussianDistribution(),
        TweedieDistribution(power=1.5), TweedieDistribution(power=4.5),
        GeneralizedHyperbolicSecand())
    for solver in ['irls', 'lbfgs', 'newton-cg', 'cd']:
        for family in families:
            glm = GeneralizedLinearRegressor(
                alpha=0, family=family, fit_intercept=False, solver=solver)
            res = glm.fit(X, y)
            assert_array_almost_equal(res.coef_, coef)


def test_glm_log_regression():
    """Test GLM regression with log link on a simple dataset
    """
    coef = [1, 2]
    X = np.array([[1, 1, 1, 1, 1], [0, 1, 2, 3, 4]]).T
    y = np.exp(np.dot(X, coef))
    families = (
        NormalDistribution(), PoissonDistribution(),
        GammaDistribution(), InverseGaussianDistribution(),
        TweedieDistribution(power=1.5), TweedieDistribution(power=4.5),
        GeneralizedHyperbolicSecand())
    for solver in ['irls', 'lbfgs', 'newton-cg']:
        for family in families:
            glm = GeneralizedLinearRegressor(
                alpha=0, family=family, link=LogLink(), fit_intercept=False,
                solver=solver, start_params='least_squares')
            res = glm.fit(X, y)
            assert_array_almost_equal(res.coef_, coef)


@pytest.mark.filterwarnings('ignore::DeprecationWarning')
def test_normal_ridge():
    """Test ridge regression for Normal distributions

    Compare to test_ridge in test_ridge.py.
    """
    rng = np.random.RandomState(0)
    alpha = 1.0

    # 1. With more samples than features
    n_samples, n_features, n_predict = 10, 5, 10
    y = rng.randn(n_samples)
    X = rng.randn(n_samples, n_features)
    T = rng.randn(n_predict, n_features)

    # GLM has 1/(2*n) * Loss + 1/2*L2, Ridge has Loss + L2
    ridge = Ridge(alpha=alpha*n_samples, fit_intercept=True, tol=1e-6,
                  solver='svd', normalize=False)
    ridge.fit(X, y)
    for solver in ['irls', 'lbfgs', 'newton-cg', 'cd']:
        glm = GeneralizedLinearRegressor(alpha=1.0, l1_ratio=0,
                                         family='normal', link='identity',
                                         fit_intercept=True, tol=1e-6,
                                         max_iter=100, solver=solver,
                                         random_state=42)
        glm.fit(X, y)
        assert_equal(glm.coef_.shape, (X.shape[1], ))
        assert_array_almost_equal(glm.coef_, ridge.coef_)
        assert_almost_equal(glm.intercept_, ridge.intercept_)
        assert_array_almost_equal(glm.predict(T), ridge.predict(T))

    ridge = Ridge(alpha=alpha*n_samples, fit_intercept=False, tol=1e-6,
                  solver='svd', normalize=False)
    ridge.fit(X, y)
    glm = GeneralizedLinearRegressor(alpha=1.0, l1_ratio=0, tol=1e-6,
                                     family='normal', link='identity',
                                     fit_intercept=False, solver='irls',
                                     fit_dispersion='chisqr')
    glm.fit(X, y)
    assert_equal(glm.coef_.shape, (X.shape[1], ))
    assert_array_almost_equal(glm.coef_, ridge.coef_)
    assert_almost_equal(glm.intercept_, ridge.intercept_)
    assert_array_almost_equal(glm.predict(T), ridge.predict(T))
    mu = glm.predict(X)
    assert_almost_equal(glm.dispersion_,
                        np.sum((y-mu)**2/(n_samples-n_features)))

    # 2. With more features than samples and sparse
    n_samples, n_features, n_predict = 5, 10, 10
    y = rng.randn(n_samples)
    X = sparse.csr_matrix(rng.randn(n_samples, n_features))
    T = sparse.csr_matrix(rng.randn(n_predict, n_features))

    # GLM has 1/(2*n) * Loss + 1/2*L2, Ridge has Loss + L2
    ridge = Ridge(alpha=alpha*n_samples, fit_intercept=True, tol=1e-9,
                  solver='sag', normalize=False, max_iter=100000)
    ridge.fit(X, y)
    for solver in ['irls', 'lbfgs', 'newton-cg', 'cd']:
        glm = GeneralizedLinearRegressor(alpha=1.0, l1_ratio=0, tol=1e-8,
                                         family='normal', link='identity',
                                         fit_intercept=True, solver=solver,
                                         max_iter=300, random_state=42)
        glm.fit(X, y)
        assert_equal(glm.coef_.shape, (X.shape[1], ))
        assert_array_almost_equal(glm.coef_, ridge.coef_, decimal=5)
        assert_almost_equal(glm.intercept_, ridge.intercept_, decimal=5)
        assert_array_almost_equal(glm.predict(T), ridge.predict(T), decimal=5)

    ridge = Ridge(alpha=alpha*n_samples, fit_intercept=False, tol=1e-7,
                  solver='sag', normalize=False, max_iter=1000)
    ridge.fit(X, y)
    glm = GeneralizedLinearRegressor(alpha=1.0, l1_ratio=0, tol=1e-7,
                                     family='normal', link='identity',
                                     fit_intercept=False, solver='irls')
    glm.fit(X, y)
    assert_equal(glm.coef_.shape, (X.shape[1], ))
    assert_array_almost_equal(glm.coef_, ridge.coef_)
    assert_almost_equal(glm.intercept_, ridge.intercept_)
    assert_array_almost_equal(glm.predict(T), ridge.predict(T))


def test_poisson_ridge():
    """Test ridge regression with poisson family and LogLink

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
    s_dec = {'irls': 7, 'lbfgs': 5, 'newton-cg': 5, 'cd': 7}
    s_tol = {'irls': 1e-8, 'lbfgs': 1e-7, 'newton-cg': 1e-7, 'cd': 1e-8}
    for solver in ['irls', 'lbfgs', 'newton-cg', 'cd']:
        glm = GeneralizedLinearRegressor(alpha=1, l1_ratio=0,
                                         fit_intercept=True, family='poisson',
                                         link='log', tol=s_tol[solver],
                                         solver=solver, max_iter=300,
                                         random_state=42)
        glm.fit(X, y)
        assert_almost_equal(glm.intercept_, -0.12889386979,
                            decimal=s_dec[solver])
        assert_array_almost_equal(glm.coef_, [0.29019207995, 0.03741173122],
                                  decimal=s_dec[solver])


def test_normal_enet():
    """Tet elastic net regression with normal/gaussian family"""
    rng = np.random.RandomState(0)
    alpha, l1_ratio = 0.3, 0.7
    n_samples, n_features = 20, 2
    X = rng.randn(n_samples, n_features).copy(order='F')
    beta = rng.randn(n_features)
    y = 2 + np.dot(X, beta) + rng.randn(n_samples)

    glm = GeneralizedLinearRegressor(alpha=alpha, l1_ratio=l1_ratio,
                                     family='normal', link='identity',
                                     fit_intercept=True, tol=1e-8,
                                     max_iter=100, selection='cyclic',
                                     solver='cd', start_params='zero',
                                     check_input=False)
    glm.fit(X, y)

    enet = ElasticNet(alpha=alpha, l1_ratio=l1_ratio, fit_intercept=True,
                      normalize=False, tol=1e-8, copy_X=True)
    enet.fit(X, y)

    assert_almost_equal(glm.intercept_, enet.intercept_, decimal=7)
    assert_array_almost_equal(glm.coef_, enet.coef_, decimal=7)


def test_poisson_enet():
    """Test elastic net regression with poisson family and LogLink

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
    glm = GeneralizedLinearRegressor(alpha=1, l1_ratio=0.5, family='poisson',
                                     link='log', solver='cd', tol=1e-7,
                                     selection='random', random_state=42)
    glm.fit(X, y)
    assert_almost_equal(glm.intercept_, glmnet_intercept, decimal=7)
    assert_array_almost_equal(glm.coef_, glmnet_coef, decimal=7)

    # same for start_params='zero' and selection='cyclic'
    # with reduced precision
    glm = GeneralizedLinearRegressor(alpha=1, l1_ratio=0.5, family='poisson',
                                     link='log', solver='cd', tol=1e-5,
                                     selection='cyclic', start_params='zero')
    glm.fit(X, y)
    assert_almost_equal(glm.intercept_, glmnet_intercept, decimal=4)
    assert_array_almost_equal(glm.coef_, glmnet_coef, decimal=4)

    # start_params='least_squares' with different alpha
    glm = GeneralizedLinearRegressor(alpha=0.005, l1_ratio=0.5,
                                     family='poisson',
                                     link='log', solver='cd', tol=1e-5,
                                     start_params='zero')
    glm.fit(X, y)
    # warm start with original alpha and use of sparse matrices
    glm.warm_start = True
    glm.alpha = 1
    X = sparse.csr_matrix(X)
    glm.fit(X, y)
    assert_almost_equal(glm.intercept_, glmnet_intercept, decimal=4)
    assert_array_almost_equal(glm.coef_, glmnet_coef, decimal=4)
