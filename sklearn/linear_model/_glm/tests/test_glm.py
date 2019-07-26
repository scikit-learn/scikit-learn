# Authors: Christian Lorentzen <lorentzen.ch@gmail.com>
#
# License: BSD 3 clause

import numpy as np
from numpy.testing import assert_allclose
import pytest

from sklearn.datasets import make_regression
from sklearn.linear_model import GeneralizedLinearRegressor
from sklearn.linear_model._glm.link import (
    IdentityLink,
    LogLink,
)
from sklearn.linear_model._glm.distribution import (
    TweedieDistribution,
    NormalDistribution, PoissonDistribution,
    GammaDistribution, InverseGaussianDistribution,
)
from sklearn.linear_model import Ridge
from sklearn.metrics import mean_absolute_error
from sklearn.exceptions import ConvergenceWarning

GLM_SOLVERS = ['lbfgs']


@pytest.fixture(scope="module")
def regression_data():
    X, y = make_regression(n_samples=107,
                           n_features=10,
                           n_informative=80, noise=0.5,
                           random_state=2)
    return X, y


def test_sample_weights_validation():
    """Test the raised errors in the validation of sample_weight."""
    # scalar value but not positive
    X = [[1]]
    y = [1]
    weights = 0
    glm = GeneralizedLinearRegressor(fit_intercept=False)

    # Positive weights are accepted
    glm.fit(X, y, sample_weight=1)

    # 2d array
    weights = [[0]]
    with pytest.raises(ValueError, match="must be 1D array or scalar"):
        glm.fit(X, y, weights)

    # 1d but wrong length
    weights = [1, 0]
    msg = r"sample_weight.shape == \(2,\), expected \(1,\)!"
    with pytest.raises(ValueError, match=msg):
        glm.fit(X, y, weights)


@pytest.mark.parametrize('f, fam',
                         [('normal', NormalDistribution()),
                          ('poisson', PoissonDistribution()),
                          ('gamma', GammaDistribution()),
                          ('inverse.gaussian', InverseGaussianDistribution())])
def test_glm_family_argument(f, fam):
    """Test GLM family argument set as string."""
    y = np.array([0.1, 0.5])  # in range of all distributions
    X = np.array([[1], [2]])
    glm = GeneralizedLinearRegressor(family=f, alpha=0).fit(X, y)
    assert isinstance(glm._family_instance, fam.__class__)

    glm = GeneralizedLinearRegressor(family='not a family',
                                     fit_intercept=False)
    with pytest.raises(ValueError, match="family must be"):
        glm.fit(X, y)


@pytest.mark.parametrize('l, link',
                         [('identity', IdentityLink()),
                          ('log', LogLink())])
def test_glm_link_argument(l, link):
    """Test GLM link argument set as string."""
    y = np.array([0.1, 0.5])  # in range of all distributions
    X = np.array([[1], [2]])
    glm = GeneralizedLinearRegressor(family='normal', link=l).fit(X, y)
    assert isinstance(glm._link_instance, link.__class__)

    glm = GeneralizedLinearRegressor(family='normal', link='not a link')
    with pytest.raises(ValueError, match="link must be"):
        glm.fit(X, y)


@pytest.mark.parametrize('alpha', ['not a number', -4.2])
def test_glm_alpha_argument(alpha):
    """Test GLM for invalid alpha argument."""
    y = np.array([1, 2])
    X = np.array([[1], [2]])
    glm = GeneralizedLinearRegressor(family='normal', alpha=alpha)
    with pytest.raises(ValueError,
                       match="Penalty term must be a non-negative"):
        glm.fit(X, y)


@pytest.mark.parametrize('fit_intercept', ['not bool', 1, 0, [True]])
def test_glm_fit_intercept_argument(fit_intercept):
    """Test GLM for invalid fit_intercept argument."""
    y = np.array([1, 2])
    X = np.array([[1], [1]])
    glm = GeneralizedLinearRegressor(fit_intercept=fit_intercept)
    with pytest.raises(ValueError, match="fit_intercept must be bool"):
        glm.fit(X, y)


@pytest.mark.parametrize('solver',
                         ['not a solver', 1, [1]])
def test_glm_solver_argument(solver):
    """Test GLM for invalid solver argument."""
    y = np.array([1, 2])
    X = np.array([[1], [2]])
    glm = GeneralizedLinearRegressor(solver=solver)
    with pytest.raises(ValueError):
        glm.fit(X, y)


@pytest.mark.parametrize('max_iter', ['not a number', 0, -1, 5.5, [1]])
def test_glm_max_iter_argument(max_iter):
    """Test GLM for invalid max_iter argument."""
    y = np.array([1, 2])
    X = np.array([[1], [2]])
    glm = GeneralizedLinearRegressor(max_iter=max_iter)
    with pytest.raises(ValueError, match="must be a positive integer"):
        glm.fit(X, y)


@pytest.mark.parametrize('tol', ['not a number', 0, -1.0, [1e-3]])
def test_glm_tol_argument(tol):
    """Test GLM for invalid tol argument."""
    y = np.array([1, 2])
    X = np.array([[1], [2]])
    glm = GeneralizedLinearRegressor(tol=tol)
    with pytest.raises(ValueError, match="stopping criteria must be positive"):
        glm.fit(X, y)


@pytest.mark.parametrize('warm_start', ['not bool', 1, 0, [True]])
def test_glm_warm_start_argument(warm_start):
    """Test GLM for invalid warm_start argument."""
    y = np.array([1, 2])
    X = np.array([[1], [1]])
    glm = GeneralizedLinearRegressor(warm_start=warm_start)
    with pytest.raises(ValueError, match="warm_start must be bool"):
        glm.fit(X, y)


@pytest.mark.parametrize('copy_X', ['not bool', 1, 0, [True]])
def test_glm_copy_X_argument(copy_X):
    """Test GLM for invalid copy_X arguments."""
    y = np.array([1, 2])
    X = np.array([[1], [1]])
    glm = GeneralizedLinearRegressor(copy_X=copy_X)
    with pytest.raises(ValueError, match="copy_X must be bool"):
        glm.fit(X, y)


@pytest.mark.parametrize('check_input', ['not bool', 1, 0, [True]])
def test_glm_check_input_argument(check_input):
    """Test GLM for invalid check_input argument."""
    y = np.array([1, 2])
    X = np.array([[1], [1]])
    glm = GeneralizedLinearRegressor(check_input=check_input)
    with pytest.raises(ValueError, match="check_input must be bool"):
        glm.fit(X, y)


@pytest.mark.parametrize('solver', GLM_SOLVERS)
def test_glm_identity_regression(solver):
    """Test GLM regression with identity link on a simple dataset."""
    coef = [1., 2.]
    X = np.array([[1, 1, 1, 1, 1], [0, 1, 2, 3, 4]]).T
    y = np.dot(X, coef)
    glm = GeneralizedLinearRegressor(alpha=0, family='normal', link='identity',
                                     fit_intercept=False, solver=solver,
                                     tol=1e-7)
    res = glm.fit(X, y)
    assert_allclose(res.coef_, coef, rtol=1e-6)


@pytest.mark.parametrize(
    'family',
    [NormalDistribution(), PoissonDistribution(),
     GammaDistribution(), InverseGaussianDistribution(),
     TweedieDistribution(power=1.5), TweedieDistribution(power=4.5)])
@pytest.mark.parametrize('solver, tol', [('lbfgs', 1e-6)])
def test_glm_log_regression(family, solver, tol):
    """Test GLM regression with log link on a simple dataset."""
    coef = [0.2, -0.1]
    X = np.array([[1, 1, 1, 1, 1], [0, 1, 2, 3, 4]]).T
    y = np.exp(np.dot(X, coef))
    glm = GeneralizedLinearRegressor(
                alpha=0, family=family, link='log', fit_intercept=False,
                solver=solver, tol=tol)
    res = glm.fit(X, y)
    assert_allclose(res.coef_, coef, rtol=5e-6)


@pytest.mark.parametrize('fit_intercept', [True, False])
def test_warm_start(fit_intercept):
    n_samples, n_features = 100, 10
    n_predict = 10
    X, y, coef = make_regression(n_samples=n_samples+n_predict,
                                 n_features=n_features,
                                 n_informative=n_features-2, noise=0.5,
                                 coef=True, random_state=42)

    glm1 = GeneralizedLinearRegressor(
        warm_start=False,
        fit_intercept=fit_intercept,
        max_iter=1000
    )
    glm1.fit(X, y)

    glm2 = GeneralizedLinearRegressor(
        warm_start=True,
        fit_intercept=fit_intercept,
        max_iter=1
    )
    glm2.fit(X, y)
    assert glm1.score(X, y) > glm2.score(X, y)
    glm2.set_params(max_iter=1000)
    glm2.fit(X, y)
    assert_allclose(glm1.coef_, glm2.coef_, rtol=1e-4, atol=1e-5)
    assert_allclose(glm1.score(X, y), glm2.score(X, y), rtol=1e-4)
    # TODO: investigate why this doesn't match
    # assert glm1.n_iter_ == glm2.n_iter_ + 2


@pytest.mark.parametrize('n_samples, n_features', [(100, 10), (10, 100)])
@pytest.mark.parametrize('fit_intercept', [True, False])
@pytest.mark.parametrize('solver', GLM_SOLVERS)
def test_normal_ridge_comparison(n_samples, n_features, fit_intercept, solver):
    """Test ridge regression for Normal distributions.

    Case n_samples >> n_features

    Compare to test_ridge in test_ridge.py.
    """
    alpha = 1.0
    n_predict = 10
    X, y, coef = make_regression(n_samples=n_samples+n_predict,
                                 n_features=n_features,
                                 n_informative=n_features-2, noise=0.5,
                                 coef=True, random_state=42)
    y = y[0:n_samples]
    X, T = X[0:n_samples], X[n_samples:]

    if n_samples > n_features:
        ridge_params = {"solver": "svd"}
    else:
        ridge_params = {"solver": "sag", "max_iter": 10000, "tol": 1e-9}

    # GLM has 1/(2*n) * Loss + 1/2*L2, Ridge has Loss + L2
    ridge = Ridge(alpha=alpha*n_samples, normalize=False,
                  random_state=42, **ridge_params)
    ridge.fit(X, y)

    glm = GeneralizedLinearRegressor(alpha=1.0, family='normal',
                                     link='identity', fit_intercept=True,
                                     max_iter=300, solver=solver, tol=1e-6,
                                     check_input=False)
    glm.fit(X, y)
    assert glm.coef_.shape == (X.shape[1], )
    assert_allclose(glm.coef_, ridge.coef_, rtol=5e-6)
    assert_allclose(glm.intercept_, ridge.intercept_, rtol=1e-5)
    assert_allclose(glm.predict(T), ridge.predict(T), rtol=1e-5)


@pytest.mark.parametrize('solver, tol', [('lbfgs', 1e-7)])
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
    glm = GeneralizedLinearRegressor(alpha=1,
                                     fit_intercept=True, family='poisson',
                                     link='log', tol=1e-7,
                                     solver=solver, max_iter=300,
                                     )
    glm.fit(X, y)
    assert_allclose(glm.intercept_, -0.12889386979, rtol=1e-5)
    assert_allclose(glm.coef_, [0.29019207995, 0.03741173122], rtol=1e-5)


@pytest.mark.parametrize(
        "params",
        [
            {"solver": "lbfgs"},
        ],
        ids=lambda params: ', '.join("%s=%s" % (key, val)
                                     for key,  val in params.items())
)
def test_solver_equivalence(params, regression_data):
    X, y = regression_data
    est_ref = GeneralizedLinearRegressor()
    est_ref.fit(X, y)

    estimator = GeneralizedLinearRegressor(**params)

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

    est1 = GeneralizedLinearRegressor()
    est1.fit(X, y)
    assert not hasattr(est1, "dispersion_")

    est2 = GeneralizedLinearRegressor(fit_dispersion="chisqr")
    est2.fit(X, y)
    assert isinstance(est2.dispersion_, float)

    est3 = GeneralizedLinearRegressor(fit_dispersion="deviance")
    est3.fit(X, y)
    assert isinstance(est3.dispersion_, float)

    assert_allclose(est2.dispersion_,  est3.dispersion_)


@pytest.mark.parametrize("solver", GLM_SOLVERS)
def test_convergence_warning(solver, regression_data):
    X, y = regression_data

    est = GeneralizedLinearRegressor(solver=solver,
                                     max_iter=1, tol=1e-20)
    with pytest.warns(ConvergenceWarning):
        est.fit(X, y)
