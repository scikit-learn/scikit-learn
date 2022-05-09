# Authors: Christian Lorentzen <lorentzen.ch@gmail.com>
#
# License: BSD 3 clause

from functools import partial
import re
import warnings

import numpy as np
from numpy.testing import assert_allclose
import pytest
from scipy import linalg
from scipy.optimize import root

from sklearn.base import clone
from sklearn._loss import HalfBinomialLoss
from sklearn._loss.glm_distribution import TweedieDistribution
from sklearn._loss.link import IdentityLink, LogLink

from sklearn.datasets import make_low_rank_matrix, make_regression
from sklearn.linear_model import (
    GammaRegressor,
    PoissonRegressor,
    Ridge,
    TweedieRegressor,
)
from sklearn.linear_model._glm import _GeneralizedLinearRegressor
from sklearn.linear_model._linear_loss import LinearModelLoss
from sklearn.exceptions import ConvergenceWarning
from sklearn.metrics import d2_tweedie_score
from sklearn.model_selection import train_test_split


SOLVERS = ["lbfgs", "newton-cholesky", "newton-qr-cholesky"]


class BinomialRegressor(_GeneralizedLinearRegressor):
    def _get_loss(self):
        return HalfBinomialLoss()


@pytest.fixture(scope="module")
def regression_data():
    X, y = make_regression(
        n_samples=107, n_features=10, n_informative=80, noise=0.5, random_state=2
    )
    return X, y


@pytest.fixture(
    params=zip(
        ["long", "wide"],
        [
            BinomialRegressor(),
            PoissonRegressor(),
            GammaRegressor(),
            TweedieRegressor(power=3.0),
            TweedieRegressor(power=0, link="log"),
            TweedieRegressor(power=1.5),
        ],
    )
)
def glm_dataset(global_random_seed, request):
    """Dataset with GLM solutions, well conditioned X.

    This is inspired by ols_ridge_dataset in test_ridge.py.

    The construction is based on the SVD decomposition of X = U S V'.

    Parameters
    ----------
    type : {"long", "wide"}
        If "long", then n_samples > n_features.
        If "wide", then n_features > n_samples.
    model : a GLM model

    For "wide", we return the minimum norm solution w = X' (XX')^-1 y:

        min ||w||_2 subject to X w = y

    Returns
    -------
    model : GLM model
    X : ndarray
        Last column of 1, i.e. intercept.
    y : ndarray
    coef_unpenalized : ndarray
        Minimum norm solutions, i.e. min sum(loss(w)) (with mininum ||w||_2 in
        case of ambiguity)
        Last coefficient is intercept.
    coef_penalized : ndarray
        GLM solution with alpha=l2_reg_strength=1, i.e.
        min 1/n * sum(loss) + ||w||_2^2.
        Last coefficient is intercept.
    """
    type, model = request.param
    # Make larger dim more than double as big as the smaller one.
    # This helps when constructing singular matrices like (X, X).
    if type == "long":
        n_samples, n_features = 12, 4
    else:
        n_samples, n_features = 4, 12
    k = min(n_samples, n_features)
    rng = np.random.RandomState(global_random_seed)
    X = make_low_rank_matrix(
        n_samples=n_samples,
        n_features=n_features,
        effective_rank=k,
        tail_strength=0.1,
        random_state=rng,
    )
    X[:, -1] = 1  # last columns acts as intercept
    U, s, Vt = linalg.svd(X)
    assert np.all(s) > 1e-3  # to be sure
    U1, _ = U[:, :k], U[:, k:]
    Vt1, _ = Vt[:k, :], Vt[k:, :]

    if request.param == "long":
        coef_unpenalized = rng.uniform(low=1, high=3, size=n_features)
        coef_unpenalized *= rng.choice([-1, 1], size=n_features)
        raw_prediction = X @ coef_unpenalized
    else:
        raw_prediction = rng.uniform(low=-3, high=3, size=n_samples)
        # w = X'(XX')^-1 y = V s^-1 U' y
        coef_unpenalized = Vt1.T @ np.diag(1 / s) @ U1.T @ raw_prediction

    linear_loss = LinearModelLoss(base_loss=model._get_loss(), fit_intercept=True)
    sw = np.full(shape=n_samples, fill_value=1 / n_samples)
    y = linear_loss.base_loss.link.inverse(raw_prediction)

    # Add penalty l2_reg_strength * ||coef||_2^2 for l2_reg_strength=1 and solve with
    # optimizer. Note that the problem is well conditioned such that we get accurate
    # results.
    l2_reg_strength = 1
    fun = partial(
        linear_loss.gradient,
        X=X[:, :-1],
        y=y,
        sample_weight=sw,
        l2_reg_strength=l2_reg_strength,
    )
    # Note: Toot finding is more precise then minimizing a function.
    res = root(
        fun,
        coef_unpenalized,
        method="lm",
        options={"ftol": 1e-14, "xtol": 1e-14, "gtol": 1e-14},
    )
    coef_penalized_with_intercept = res.x

    linear_loss = LinearModelLoss(base_loss=model._get_loss(), fit_intercept=False)
    fun = partial(
        linear_loss.gradient,
        X=X[:, :-1],
        y=y,
        sample_weight=sw,
        l2_reg_strength=l2_reg_strength,
    )
    res = root(
        fun,
        coef_unpenalized[:-1],
        method="lm",
        options={"ftol": 1e-14, "xtol": 1e-14, "gtol": 1e-14},
    )
    coef_penalized_without_intercept = res.x

    # To be sure
    assert np.linalg.norm(coef_penalized_with_intercept) < np.linalg.norm(
        coef_unpenalized
    )

    return (
        model,
        X,
        y,
        coef_unpenalized,
        coef_penalized_with_intercept,
        coef_penalized_without_intercept,
    )


@pytest.mark.parametrize("solver", SOLVERS)
@pytest.mark.parametrize("fit_intercept", [False, True])
def test_glm_regression(solver, fit_intercept, glm_dataset):
    """Test that GLM converges for all solvers to correct solution.

    We work with a simple constructed data set with known solution.
    """
    model, X, y, _, coef_with_intercept, coef_without_intercept = glm_dataset
    alpha = 1.0  # because glm_dataset uses this.
    params = dict(
        alpha=alpha,
        fit_intercept=fit_intercept,
        solver=solver,
        tol=1e-12,
        max_iter=1000,
    )

    model = clone(model).set_params(**params)
    X = X[:, :-1]  # remove intercept
    if fit_intercept:
        coef = coef_with_intercept
        intercept = coef[-1]
        coef = coef[:-1]
    else:
        coef = coef_without_intercept
        intercept = 0
    model.fit(X, y)

    rtol = 3e-5 if solver == "lbfgs" else 1e-11
    assert model.intercept_ == pytest.approx(intercept, rel=rtol)
    assert_allclose(model.coef_, coef, rtol=rtol)

    # Same with sample_weight.
    model = (
        clone(model).set_params(**params).fit(X, y, sample_weight=np.ones(X.shape[0]))
    )
    assert model.intercept_ == pytest.approx(intercept, rel=rtol)
    assert_allclose(model.coef_, coef, rtol=rtol)


@pytest.mark.parametrize("solver", SOLVERS)
@pytest.mark.parametrize("fit_intercept", [True, False])
def test_glm_regression_hstacked_X(solver, fit_intercept, glm_dataset):
    """Test that GLM converges for all solvers to correct solution on hstacked data.

    We work with a simple constructed data set with known solution.
    Fit on [X] with alpha is the same as fit on [X, X]/2 with alpha/2.
    For long X, [X, X] is a singular matrix.
    """
    model, X, y, _, coef_with_intercept, coef_without_intercept = glm_dataset
    n_samples, n_features = X.shape
    alpha = 1.0  # because glm_dataset uses this.
    params = dict(
        alpha=alpha / 2,
        fit_intercept=fit_intercept,
        solver=solver,
        tol=1e-12,
        max_iter=1000,
    )

    if (
        solver == "lbfgs"
        and fit_intercept is False
        and (
            isinstance(model, BinomialRegressor)
            or (isinstance(model, PoissonRegressor) and n_features > n_samples)
        )
    ):
        # Line search cannot locate an adequate point after MAXLS
        #  function and gradient evaluations.
        #  Previous x, f and g restored.
        # Possible causes: 1 error in function or gradient evaluation;
        #                  2 rounding error dominate computation.
        pytest.xfail()

    model = clone(model).set_params(**params)
    X = X[:, :-1]  # remove intercept
    X = 0.5 * np.concatenate((X, X), axis=1)
    assert np.linalg.matrix_rank(X) <= min(n_samples, n_features - 1)
    if fit_intercept:
        coef = coef_with_intercept
        intercept = coef[-1]
        coef = coef[:-1]
    else:
        coef = coef_without_intercept
        intercept = 0
    model.fit(X, y)

    rtol = 3e-5 if solver == "lbfgs" else 1e-11
    assert model.intercept_ == pytest.approx(intercept, rel=rtol)
    assert_allclose(model.coef_, np.r_[coef, coef], rtol=rtol)


@pytest.mark.parametrize("solver", SOLVERS)
@pytest.mark.parametrize("fit_intercept", [True, False])
def test_glm_regression_vstacked_X(solver, fit_intercept, glm_dataset):
    """Test that GLM converges for all solvers to correct solution on vstacked data.

    We work with a simple constructed data set with known solution.
    Fit on [X] with alpha is the same as fit on [X], [y]
                                                [X], [y] with 1 * alpha.
    It is the same alpha as the average loss stays the same.
    For wide X, [X', X'] is a singular matrix.
    """
    model, X, y, _, coef_with_intercept, coef_without_intercept = glm_dataset
    n_samples, n_features = X.shape
    alpha = 1.0  # because glm_dataset uses this.
    params = dict(
        alpha=alpha,
        fit_intercept=fit_intercept,
        solver=solver,
        tol=1e-12,
        max_iter=1000,
    )

    model = clone(model).set_params(**params)
    X = X[:, :-1]  # remove intercept
    X = np.concatenate((X, X), axis=0)
    assert np.linalg.matrix_rank(X) <= min(n_samples, n_features)
    y = np.r_[y, y]
    if fit_intercept:
        coef = coef_with_intercept
        intercept = coef[-1]
        coef = coef[:-1]
    else:
        coef = coef_without_intercept
        intercept = 0
    model.fit(X, y)

    rtol = 3e-5 if solver == "lbfgs" else 1e-11
    assert model.intercept_ == pytest.approx(intercept, rel=rtol)
    assert_allclose(model.coef_, coef, rtol=rtol)


def test_sample_weights_validation():
    """Test the raised errors in the validation of sample_weight."""
    # scalar value but not positive
    X = [[1]]
    y = [1]
    weights = 0
    glm = _GeneralizedLinearRegressor()

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


@pytest.mark.parametrize("fit_intercept", ["not bool", 1, 0, [True]])
def test_glm_fit_intercept_argument(fit_intercept):
    """Test GLM for invalid fit_intercept argument."""
    y = np.array([1, 2])
    X = np.array([[1], [1]])
    glm = _GeneralizedLinearRegressor(fit_intercept=fit_intercept)
    with pytest.raises(ValueError, match="fit_intercept must be bool"):
        glm.fit(X, y)


@pytest.mark.parametrize("solver", ["not a solver", 1, [1]])
def test_glm_solver_argument(solver):
    """Test GLM for invalid solver argument."""
    y = np.array([1, 2])
    X = np.array([[1], [2]])
    glm = _GeneralizedLinearRegressor(solver=solver)
    with pytest.raises(ValueError):
        glm.fit(X, y)


@pytest.mark.parametrize(
    "Estimator",
    [_GeneralizedLinearRegressor, PoissonRegressor, GammaRegressor, TweedieRegressor],
)
@pytest.mark.parametrize(
    "params, err_type, err_msg",
    [
        ({"max_iter": 0}, ValueError, "max_iter == 0, must be >= 1"),
        ({"max_iter": -1}, ValueError, "max_iter == -1, must be >= 1"),
        (
            {"max_iter": "not a number"},
            TypeError,
            "max_iter must be an instance of int, not str",
        ),
        (
            {"max_iter": [1]},
            TypeError,
            "max_iter must be an instance of int, not list",
        ),
        (
            {"max_iter": 5.5},
            TypeError,
            "max_iter must be an instance of int, not float",
        ),
        ({"alpha": -1}, ValueError, "alpha == -1, must be >= 0.0"),
        (
            {"alpha": "1"},
            TypeError,
            "alpha must be an instance of float, not str",
        ),
        ({"tol": -1.0}, ValueError, "tol == -1.0, must be > 0."),
        ({"tol": 0.0}, ValueError, "tol == 0.0, must be > 0.0"),
        ({"tol": 0}, ValueError, "tol == 0, must be > 0.0"),
        (
            {"tol": "1"},
            TypeError,
            "tol must be an instance of float, not str",
        ),
        (
            {"tol": [1e-3]},
            TypeError,
            "tol must be an instance of float, not list",
        ),
        ({"verbose": -1}, ValueError, "verbose == -1, must be >= 0."),
        (
            {"verbose": "1"},
            TypeError,
            "verbose must be an instance of int, not str",
        ),
        (
            {"verbose": 1.0},
            TypeError,
            "verbose must be an instance of int, not float",
        ),
    ],
)
def test_glm_scalar_argument(Estimator, params, err_type, err_msg):
    """Test GLM for invalid parameter arguments."""
    y = np.array([1, 2])
    X = np.array([[1], [2]])
    glm = Estimator(**params)
    with pytest.raises(err_type, match=err_msg):
        glm.fit(X, y)


@pytest.mark.parametrize("warm_start", ["not bool", 1, 0, [True]])
def test_glm_warm_start_argument(warm_start):
    """Test GLM for invalid warm_start argument."""
    y = np.array([1, 2])
    X = np.array([[1], [1]])
    glm = _GeneralizedLinearRegressor(warm_start=warm_start)
    with pytest.raises(ValueError, match="warm_start must be bool"):
        glm.fit(X, y)


@pytest.mark.parametrize(
    "glm",
    [
        TweedieRegressor(power=3),
        PoissonRegressor(),
        GammaRegressor(),
        TweedieRegressor(power=1.5),
    ],
)
def test_glm_wrong_y_range(glm):
    y = np.array([-1, 2])
    X = np.array([[1], [1]])
    msg = r"Some value\(s\) of y are out of the valid range of the loss"
    with pytest.raises(ValueError, match=msg):
        glm.fit(X, y)


@pytest.mark.parametrize("fit_intercept", [False, True])
def test_glm_identity_regression(fit_intercept):
    """Test GLM regression with identity link on a simple dataset."""
    coef = [1.0, 2.0]
    X = np.array([[1, 1, 1, 1, 1], [0, 1, 2, 3, 4]]).T
    y = np.dot(X, coef)
    glm = _GeneralizedLinearRegressor(
        alpha=0,
        fit_intercept=fit_intercept,
        tol=1e-12,
    )
    if fit_intercept:
        glm.fit(X[:, 1:], y)
        assert_allclose(glm.coef_, coef[1:], rtol=1e-10)
        assert_allclose(glm.intercept_, coef[0], rtol=1e-10)
    else:
        glm.fit(X, y)
        assert_allclose(glm.coef_, coef, rtol=1e-12)


@pytest.mark.parametrize("fit_intercept", [False, True])
@pytest.mark.parametrize("alpha", [0.0, 1.0])
@pytest.mark.parametrize(
    "GLMEstimator", [_GeneralizedLinearRegressor, PoissonRegressor, GammaRegressor]
)
def test_glm_sample_weight_consistency(fit_intercept, alpha, GLMEstimator):
    """Test that the impact of sample_weight is consistent"""
    rng = np.random.RandomState(0)
    n_samples, n_features = 10, 5

    X = rng.rand(n_samples, n_features)
    y = rng.rand(n_samples)
    glm_params = dict(alpha=alpha, fit_intercept=fit_intercept)

    glm = GLMEstimator(**glm_params).fit(X, y)
    coef = glm.coef_.copy()

    # sample_weight=np.ones(..) should be equivalent to sample_weight=None
    sample_weight = np.ones(y.shape)
    glm.fit(X, y, sample_weight=sample_weight)
    assert_allclose(glm.coef_, coef, rtol=1e-12)

    # sample_weight are normalized to 1 so, scaling them has no effect
    sample_weight = 2 * np.ones(y.shape)
    glm.fit(X, y, sample_weight=sample_weight)
    assert_allclose(glm.coef_, coef, rtol=1e-12)

    # setting one element of sample_weight to 0 is equivalent to removing
    # the corresponding sample
    sample_weight = np.ones(y.shape)
    sample_weight[-1] = 0
    glm.fit(X, y, sample_weight=sample_weight)
    coef1 = glm.coef_.copy()
    glm.fit(X[:-1], y[:-1])
    assert_allclose(glm.coef_, coef1, rtol=1e-12)

    # check that multiplying sample_weight by 2 is equivalent
    # to repeating corresponding samples twice
    X2 = np.concatenate([X, X[: n_samples // 2]], axis=0)
    y2 = np.concatenate([y, y[: n_samples // 2]])
    sample_weight_1 = np.ones(len(y))
    sample_weight_1[: n_samples // 2] = 2

    glm1 = GLMEstimator(**glm_params).fit(X, y, sample_weight=sample_weight_1)

    glm2 = GLMEstimator(**glm_params).fit(X2, y2, sample_weight=None)
    assert_allclose(glm1.coef_, glm2.coef_)


@pytest.mark.parametrize("solver", SOLVERS)
@pytest.mark.parametrize("fit_intercept", [True, False])
@pytest.mark.parametrize(
    "estimator",
    [
        PoissonRegressor(),
        GammaRegressor(),
        TweedieRegressor(power=3.0),
        TweedieRegressor(power=0, link="log"),
        TweedieRegressor(power=1.5),
        TweedieRegressor(power=4.5),
    ],
)
def test_glm_log_regression(solver, fit_intercept, estimator):
    """Test GLM regression with log link on a simple dataset."""
    coef = [0.2, -0.1]
    X = np.array([[0, 1, 2, 3, 4], [1, 1, 1, 1, 1]]).T
    y = np.exp(np.dot(X, coef))
    glm = clone(estimator).set_params(
        alpha=0,
        fit_intercept=fit_intercept,
        solver=solver,
        tol=1e-8,
    )
    if fit_intercept:
        res = glm.fit(X[:, :-1], y)
        assert_allclose(res.coef_, coef[:-1], rtol=1e-6)
        assert_allclose(res.intercept_, coef[-1], rtol=1e-6)
    else:
        res = glm.fit(X, y)
        assert_allclose(res.coef_, coef, rtol=2e-6)


@pytest.mark.parametrize("solver", SOLVERS)
@pytest.mark.parametrize("fit_intercept", [True, False])
def test_warm_start(solver, fit_intercept, global_random_seed):
    n_samples, n_features = 100, 10
    X, y = make_regression(
        n_samples=n_samples,
        n_features=n_features,
        n_informative=n_features - 2,
        bias=fit_intercept * 1.0,
        noise=1.0,
        random_state=global_random_seed,
    )
    y = np.abs(y)  # Poisson requires non-negative targets.
    params = {"solver": solver, "fit_intercept": fit_intercept, "tol": 1e-10}

    glm1 = PoissonRegressor(warm_start=False, max_iter=1000, **params)
    glm1.fit(X, y)

    glm2 = PoissonRegressor(warm_start=True, max_iter=1, **params)
    # As we intentionally set max_iter=1, the solver will issue a
    # ConvergenceWarning which we here simply ignore.
    with warnings.catch_warnings():
        warnings.filterwarnings("ignore", category=ConvergenceWarning)
        glm2.fit(X, y)

    linear_loss = LinearModelLoss(
        base_loss=glm1._get_loss(),
        fit_intercept=fit_intercept,
    )
    sw = np.full_like(y, fill_value=1 / n_samples)

    objective_glm1 = linear_loss.loss(
        coef=np.r_[glm1.coef_, glm1.intercept_] if fit_intercept else glm1.coef_,
        X=X,
        y=y,
        sample_weight=sw,
        l2_reg_strength=1.0,
    )
    objective_glm2 = linear_loss.loss(
        coef=np.r_[glm2.coef_, glm2.intercept_] if fit_intercept else glm2.coef_,
        X=X,
        y=y,
        sample_weight=sw,
        l2_reg_strength=1.0,
    )
    assert objective_glm1 < objective_glm2

    glm2.set_params(max_iter=1000)
    glm2.fit(X, y)
    # The two models are not exactly identical since the lbfgs solver
    # computes the approximate hessian from previous iterations, which
    # will not be strictly identical in the case of a warm start.
    rtol = 2e-4 if solver == "lbfgs" else 1e-9
    assert_allclose(glm1.coef_, glm2.coef_, rtol=rtol)
    assert_allclose(glm1.score(X, y), glm2.score(X, y), rtol=1e-5)


# FIXME: 'normalize' to be removed in 1.2 in LinearRegression
@pytest.mark.filterwarnings("ignore:'normalize' was deprecated")
@pytest.mark.parametrize("n_samples, n_features", [(100, 10), (10, 100)])
@pytest.mark.parametrize("fit_intercept", [True, False])
@pytest.mark.parametrize("sample_weight", [None, True])
def test_normal_ridge_comparison(
    n_samples, n_features, fit_intercept, sample_weight, request
):
    """Compare with Ridge regression for Normal distributions."""
    test_size = 10
    X, y = make_regression(
        n_samples=n_samples + test_size,
        n_features=n_features,
        n_informative=n_features - 2,
        noise=0.5,
        random_state=42,
    )

    if n_samples > n_features:
        ridge_params = {"solver": "svd"}
    else:
        ridge_params = {"solver": "saga", "max_iter": 1000000, "tol": 1e-7}

    (
        X_train,
        X_test,
        y_train,
        y_test,
    ) = train_test_split(X, y, test_size=test_size, random_state=0)

    alpha = 1.0
    if sample_weight is None:
        sw_train = None
        alpha_ridge = alpha * n_samples
    else:
        sw_train = np.random.RandomState(0).rand(len(y_train))
        alpha_ridge = alpha * sw_train.sum()

    # GLM has 1/(2*n) * Loss + 1/2*L2, Ridge has Loss + L2
    ridge = Ridge(
        alpha=alpha_ridge,
        normalize=False,
        random_state=42,
        fit_intercept=fit_intercept,
        **ridge_params,
    )
    ridge.fit(X_train, y_train, sample_weight=sw_train)

    glm = _GeneralizedLinearRegressor(
        alpha=alpha,
        fit_intercept=fit_intercept,
        max_iter=300,
        tol=1e-5,
    )
    glm.fit(X_train, y_train, sample_weight=sw_train)
    assert glm.coef_.shape == (X.shape[1],)
    assert_allclose(glm.coef_, ridge.coef_, atol=5e-5)
    assert_allclose(glm.intercept_, ridge.intercept_, rtol=1e-5)
    assert_allclose(glm.predict(X_train), ridge.predict(X_train), rtol=2e-4)
    assert_allclose(glm.predict(X_test), ridge.predict(X_test), rtol=2e-4)


@pytest.mark.parametrize("solver", ["lbfgs", "newton-cholesky", "newton-qr-cholesky"])
def test_poisson_glmnet(solver):
    """Compare Poisson regression with L2 regularization and LogLink to glmnet"""
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
    glm = PoissonRegressor(
        alpha=1,
        fit_intercept=True,
        tol=1e-7,
        max_iter=300,
        solver=solver,
    )
    glm.fit(X, y)
    assert_allclose(glm.intercept_, -0.12889386979, rtol=1e-5)
    assert_allclose(glm.coef_, [0.29019207995, 0.03741173122], rtol=1e-5)


def test_convergence_warning(regression_data):
    X, y = regression_data

    est = _GeneralizedLinearRegressor(max_iter=1, tol=1e-20)
    with pytest.warns(ConvergenceWarning):
        est.fit(X, y)


@pytest.mark.parametrize(
    "name, link_class", [("identity", IdentityLink), ("log", LogLink)]
)
def test_tweedie_link_argument(name, link_class):
    """Test GLM link argument set as string."""
    y = np.array([0.1, 0.5])  # in range of all distributions
    X = np.array([[1], [2]])
    glm = TweedieRegressor(power=1, link=name).fit(X, y)
    assert isinstance(glm._base_loss.link, link_class)

    glm = TweedieRegressor(power=1, link="not a link")
    with pytest.raises(
        ValueError,
        match=re.escape("The link must be an element of ['auto', 'identity', 'log']"),
    ):
        glm.fit(X, y)


@pytest.mark.parametrize(
    "power, expected_link_class",
    [
        (0, IdentityLink),  # normal
        (1, LogLink),  # poisson
        (2, LogLink),  # gamma
        (3, LogLink),  # inverse-gaussian
    ],
)
def test_tweedie_link_auto(power, expected_link_class):
    """Test that link='auto' delivers the expected link function"""
    y = np.array([0.1, 0.5])  # in range of all distributions
    X = np.array([[1], [2]])
    glm = TweedieRegressor(link="auto", power=power).fit(X, y)
    assert isinstance(glm._base_loss.link, expected_link_class)


@pytest.mark.parametrize("power", [0, 1, 1.5, 2, 3])
@pytest.mark.parametrize("link", ["log", "identity"])
def test_tweedie_score(regression_data, power, link):
    """Test that GLM score equals d2_tweedie_score for Tweedie losses."""
    X, y = regression_data
    # make y positive
    y = np.abs(y) + 1.0
    glm = TweedieRegressor(power=power, link=link).fit(X, y)
    assert glm.score(X, y) == pytest.approx(
        d2_tweedie_score(y, glm.predict(X), power=power)
    )


@pytest.mark.parametrize(
    "estimator, value",
    [
        (PoissonRegressor(), True),
        (GammaRegressor(), True),
        (TweedieRegressor(power=1.5), True),
        (TweedieRegressor(power=0), False),
    ],
)
def test_tags(estimator, value):
    assert estimator._get_tags()["requires_positive_y"] is value


# TODO(1.3): remove
@pytest.mark.parametrize(
    "est, family",
    [
        (PoissonRegressor(), "poisson"),
        (GammaRegressor(), "gamma"),
        (TweedieRegressor(), TweedieDistribution()),
        (TweedieRegressor(power=2), TweedieDistribution(power=2)),
        (TweedieRegressor(power=3), TweedieDistribution(power=3)),
    ],
)
def test_family_deprecation(est, family):
    """Test backward compatibility of the family property."""
    with pytest.warns(FutureWarning, match="`family` was deprecated"):
        if isinstance(family, str):
            assert est.family == family
        else:
            assert est.family.__class__ == family.__class__
            assert est.family.power == family.power
