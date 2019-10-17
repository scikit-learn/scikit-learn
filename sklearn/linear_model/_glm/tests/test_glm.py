# Authors: Christian Lorentzen <lorentzen.ch@gmail.com>
#
# License: BSD 3 clause

import numpy as np
from numpy.testing import assert_allclose
import pytest

from sklearn.datasets import make_regression
from sklearn.linear_model._glm import GeneralizedLinearRegressor
from sklearn.linear_model import (
    TweedieRegressor,
    PoissonRegressor,
    GammaRegressor
)
from sklearn.linear_model._glm.link import (
    IdentityLink,
    LogLink,
)
from sklearn._loss.glm_distribution import (
    TweedieDistribution,
    NormalDistribution, PoissonDistribution,
    GammaDistribution, InverseGaussianDistribution,
)
from sklearn.linear_model import Ridge
from sklearn.exceptions import ConvergenceWarning
from sklearn.model_selection import train_test_split


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


@pytest.mark.parametrize('name, instance',
                         [('normal', NormalDistribution()),
                          ('poisson', PoissonDistribution()),
                          ('gamma', GammaDistribution()),
                          ('inverse-gaussian', InverseGaussianDistribution())])
def test_glm_family_argument(name, instance):
    """Test GLM family argument set as string."""
    y = np.array([0.1, 0.5])  # in range of all distributions
    X = np.array([[1], [2]])
    glm = GeneralizedLinearRegressor(family=name, alpha=0).fit(X, y)
    assert isinstance(glm._family_instance, instance.__class__)

    glm = GeneralizedLinearRegressor(family='not a family',
                                     fit_intercept=False)
    with pytest.raises(ValueError, match="family must be"):
        glm.fit(X, y)


@pytest.mark.parametrize('name, instance',
                         [('identity', IdentityLink()),
                          ('log', LogLink())])
def test_glm_link_argument(name, instance):
    """Test GLM link argument set as string."""
    y = np.array([0.1, 0.5])  # in range of all distributions
    X = np.array([[1], [2]])
    glm = GeneralizedLinearRegressor(family='normal', link=name).fit(X, y)
    assert isinstance(glm._link_instance, instance.__class__)

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


def test_glm_identity_regression():
    """Test GLM regression with identity link on a simple dataset."""
    coef = [1., 2.]
    X = np.array([[1, 1, 1, 1, 1], [0, 1, 2, 3, 4]]).T
    y = np.dot(X, coef)
    glm = GeneralizedLinearRegressor(alpha=0, family='normal', link='identity',
                                     fit_intercept=False)
    glm.fit(X, y)
    assert_allclose(glm.coef_, coef, rtol=1e-6)


def test_glm_sample_weight_consistentcy():
    """Test that the impact of sample_weight is consistent"""
    rng = np.random.RandomState(0)
    n_samples, n_features = 10, 5

    X = rng.rand(n_samples, n_features)
    y = rng.rand(n_samples)
    glm = GeneralizedLinearRegressor(alpha=0, family='normal', link='identity',
                                     fit_intercept=False)
    glm.fit(X, y)
    coef = glm.coef_.copy()

    # sample_weight=np.ones(..) should be equivalent to sample_weight=None
    sample_weight = np.ones(y.shape)
    glm.fit(X, y, sample_weight=sample_weight)
    assert_allclose(glm.coef_, coef, rtol=1e-6)

    # sample_weight are normalized to 1 so, scaling them has no effect
    sample_weight = 2*np.ones(y.shape)
    glm.fit(X, y, sample_weight=sample_weight)
    assert_allclose(glm.coef_, coef, rtol=1e-6)

    # setting one element of sample_weight to 0 is equivalent to removing
    # the correspoding sample
    sample_weight = np.ones(y.shape)
    sample_weight[-1] = 0
    glm.fit(X, y, sample_weight=sample_weight)
    coef1 = glm.coef_.copy()
    glm.fit(X[:-1], y[:-1])
    assert_allclose(glm.coef_, coef1, rtol=1e-6)


@pytest.mark.parametrize(
    'family',
    [NormalDistribution(), PoissonDistribution(),
     GammaDistribution(), InverseGaussianDistribution(),
     TweedieDistribution(power=1.5), TweedieDistribution(power=4.5)])
def test_glm_log_regression(family):
    """Test GLM regression with log link on a simple dataset."""
    coef = [0.2, -0.1]
    X = np.array([[1, 1, 1, 1, 1], [0, 1, 2, 3, 4]]).T
    y = np.exp(np.dot(X, coef))
    glm = GeneralizedLinearRegressor(
                alpha=0, family=family, link='log', fit_intercept=False,
                tol=1e-6)
    res = glm.fit(X, y)
    assert_allclose(res.coef_, coef, rtol=5e-6)


@pytest.mark.parametrize('fit_intercept', [True, False])
def test_warm_start(fit_intercept):
    n_samples, n_features = 110, 10
    X, y, coef = make_regression(n_samples=n_samples,
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
    # The two model are not exactly identical since the lbfgs solver
    # computes the approximate hessian from previous iterations, which
    # will not be strictly identical in the case of a warm start.
    assert_allclose(glm1.coef_, glm2.coef_, rtol=1e-5)
    assert_allclose(glm1.score(X, y), glm2.score(X, y), rtol=1e-4)


@pytest.mark.parametrize('n_samples, n_features', [(100, 10), (10, 100)])
@pytest.mark.parametrize('fit_intercept', [True, False])
@pytest.mark.parametrize('sample_weight', [None, pytest.mark.xfail('rand')])
def test_normal_ridge_comparison(n_samples, n_features, fit_intercept,
                                 sample_weight, request):
    """Compare with Ridge regression for Normal distributions."""
    test_size = 10
    X, y = make_regression(n_samples=n_samples + test_size,
                           n_features=n_features,
                           n_informative=n_features-2, noise=0.5,
                           random_state=42)

    if n_samples > n_features:
        ridge_params = {"solver": "svd"}
    else:
        ridge_params = {"solver": "saga", "max_iter": 1000000, "tol": 1e-9}

    X_train, X_test, y_train, y_test, = train_test_split(
        X, y, test_size=test_size, random_state=0
    )

    if sample_weight is None:
        alpha = 1.0
        sw_train = None
    else:
        sw_train = np.random.RandomState(0).rand(len(y_train))
        alpha = 0.0
        sw_train /= sw_train.sum()
        request.applymarker(pytest.mark.xfail(
            run=False, reason=('TODO: GLM / Ridge comparison with '
                               'sample_weight should be fixed')))

    # GLM has 1/(2*n) * Loss + 1/2*L2, Ridge has Loss + L2
    ridge = Ridge(alpha=alpha*n_samples, normalize=False,
                  random_state=42, fit_intercept=fit_intercept,
                  **ridge_params)
    ridge.fit(X_train, y_train, sample_weight=sw_train)

    glm = GeneralizedLinearRegressor(alpha=1.0, family='normal',
                                     link='identity',
                                     fit_intercept=fit_intercept,
                                     max_iter=300)
    glm.fit(X_train, y_train, sample_weight=sw_train)
    assert glm.coef_.shape == (X.shape[1], )
    assert_allclose(glm.coef_, ridge.coef_, atol=5e-5)
    assert_allclose(glm.intercept_, ridge.intercept_, rtol=1e-5)
    assert_allclose(glm.predict(X_train), ridge.predict(X_train), rtol=2e-4)
    assert_allclose(glm.predict(X_test), ridge.predict(X_test), rtol=2e-4)


def test_poisson_glmnet():
    """Compare Poisson regression with L2 regularization and LogLink to glmnet
    """
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
                                     max_iter=300)
    glm.fit(X, y)
    assert_allclose(glm.intercept_, -0.12889386979, rtol=1e-5)
    assert_allclose(glm.coef_, [0.29019207995, 0.03741173122], rtol=1e-5)


def test_convergence_warning(regression_data):
    X, y = regression_data

    est = GeneralizedLinearRegressor(max_iter=1, tol=1e-20)
    with pytest.warns(ConvergenceWarning):
        est.fit(X, y)


def test_poisson_regression_family(regression_data):
    est = PoissonRegressor()
    est.family == "poisson"

    msg = "PoissonRegressor.family must be 'poisson'!"
    with pytest.raises(ValueError, match=msg):
        est.family = 0


def test_gamma_regression_family(regression_data):
    est = GammaRegressor()
    est.family == "gamma"

    msg = "GammaRegressor.family must be 'gamma'!"
    with pytest.raises(ValueError, match=msg):
        est.family = 0


def test_tweedie_regression_family(regression_data):
    power = 2.0
    est = TweedieRegressor(power=power)
    assert isinstance(est.family, TweedieDistribution)
    assert est.family.power == power
    msg = "TweedieRegressor.family must be of type TweedieDistribution!"
    with pytest.raises(TypeError, match=msg):
        est.family = None


@pytest.mark.parametrize(
        'estimator, value',
        [
            (PoissonRegressor(), True),
            (GammaRegressor(), True),
            (TweedieRegressor(power=1.5), True),
            (TweedieRegressor(power=0), False)
        ],
)
def test_tags(estimator, value):
    assert estimator._get_tags()['requires_positive_y'] is value
