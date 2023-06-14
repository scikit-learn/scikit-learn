# Authors: Christian Lorentzen <lorentzen.ch@gmail.com>
#
# License: BSD 3 clause

from functools import partial
import itertools
import warnings

import numpy as np
from numpy.testing import assert_allclose
import pytest
import scipy
from scipy import linalg
from scipy.optimize import minimize, root
from scipy.sparse.linalg import lsmr

from sklearn.base import clone
from sklearn._loss import (
    HalfBinomialLoss,
    HalfPoissonLoss,
    HalfTweedieLoss,
    HalfMultinomialLoss,
)
from sklearn._loss.link import IdentityLink, LogLink

from sklearn.datasets import make_low_rank_matrix, make_regression, make_classification
from sklearn.linear_model import (
    GammaRegressor,
    PoissonRegressor,
    Ridge,
    TweedieRegressor,
)
from sklearn.linear_model._glm import _GeneralizedLinearRegressor
from sklearn.linear_model._glm._newton_solver import (
    NewtonCholeskySolver,
    NewtonLSMRSolver,
)
from sklearn.linear_model._linear_loss import (
    LinearModelLoss,
    Multinomial_LDL_Decomposition,
)
from sklearn.exceptions import ConvergenceWarning
from sklearn.metrics import d2_tweedie_score, mean_poisson_deviance
from sklearn.model_selection import train_test_split
from sklearn.utils.fixes import parse_version, sp_version


SOLVERS = ["lbfgs", "newton-cholesky", "newton-lsmr"]


class BinomialRegressor(_GeneralizedLinearRegressor):
    def _get_loss(self):
        return HalfBinomialLoss()


def _special_minimize(fun, grad, x, tol_NM, tol):
    # Find good starting point by Nelder-Mead
    res_NM = minimize(
        fun, x, method="Nelder-Mead", options={"xatol": tol_NM, "fatol": tol_NM}
    )
    # Now refine via root finding on the gradient of the function, which is
    # more precise than minimizing the function itself.
    res = root(
        grad,
        res_NM.x,
        method="lm",
        options={"ftol": tol, "xtol": tol, "gtol": tol},
    )
    return res.x


@pytest.fixture(scope="module")
def regression_data():
    X, y = make_regression(
        n_samples=107, n_features=10, n_informative=80, noise=0.5, random_state=2
    )
    return X, y


@pytest.fixture(
    params=itertools.product(
        ["long", "wide"],
        [
            BinomialRegressor(),
            PoissonRegressor(),
            GammaRegressor(),
            # TweedieRegressor(power=3.0),  # too difficult
            # TweedieRegressor(power=0, link="log"),  # too difficult
            TweedieRegressor(power=1.5),
        ],
    ),
    ids=lambda param: f"{param[0]}-{param[1]}",
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

    For "wide", we return the minimum norm solution:

        min ||w||_2 subject to w = argmin deviance(X, y, w)

    Note that the deviance is always minimized if y = inverse_link(X w) is possible to
    achieve, which it is in the wide data case. Therefore, we can construct the
    solution with minimum norm like (wide) OLS:

        min ||w||_2 subject to link(y) = raw_prediction = X w

    Returns
    -------
    model : GLM model
    X : ndarray
        Last column of 1, i.e. intercept.
    y : ndarray
    coef_unpenalized : ndarray
        Minimum norm solutions, i.e. min sum(loss(w)) (with minimum ||w||_2 in
        case of ambiguity)
        Last coefficient is intercept.
    coef_penalized : ndarray
        GLM solution with alpha=l2_reg_strength=1, i.e.
        min 1/n * sum(loss) + ||w[:-1]||_2^2.
        Last coefficient is intercept.
    l2_reg_strength : float
        Always equal 1.
    """
    data_type, model = request.param
    # Make larger dim more than double as big as the smaller one.
    # This helps when constructing singular matrices like (X, X).
    if data_type == "long":
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
    U, s, Vt = linalg.svd(X, full_matrices=False)
    assert np.all(s > 1e-3)  # to be sure
    assert np.max(s) / np.min(s) < 100  # condition number of X

    if data_type == "long":
        coef_unpenalized = rng.uniform(low=1, high=3, size=n_features)
        coef_unpenalized *= rng.choice([-1, 1], size=n_features)
        raw_prediction = X @ coef_unpenalized
    else:
        raw_prediction = rng.uniform(low=-3, high=3, size=n_samples)
        # minimum norm solution min ||w||_2 such that raw_prediction = X w:
        # w = X'(XX')^-1 raw_prediction = V s^-1 U' raw_prediction
        coef_unpenalized = Vt.T @ np.diag(1 / s) @ U.T @ raw_prediction

    linear_loss = LinearModelLoss(base_loss=model._get_loss(), fit_intercept=True)
    sw = np.full(shape=n_samples, fill_value=1 / n_samples)
    y = linear_loss.base_loss.link.inverse(raw_prediction)

    # Add penalty l2_reg_strength * ||coef||_2^2 for l2_reg_strength=1 and solve with
    # optimizer. Note that the problem is well conditioned such that we get accurate
    # results.
    l2_reg_strength = 1
    fun = partial(
        linear_loss.loss,
        X=X[:, :-1],
        y=y,
        sample_weight=sw,
        l2_reg_strength=l2_reg_strength,
    )
    grad = partial(
        linear_loss.gradient,
        X=X[:, :-1],
        y=y,
        sample_weight=sw,
        l2_reg_strength=l2_reg_strength,
    )
    coef_penalized_with_intercept = _special_minimize(
        fun, grad, coef_unpenalized, tol_NM=1e-6, tol=1e-14
    )

    linear_loss = LinearModelLoss(base_loss=model._get_loss(), fit_intercept=False)
    fun = partial(
        linear_loss.loss,
        X=X[:, :-1],
        y=y,
        sample_weight=sw,
        l2_reg_strength=l2_reg_strength,
    )
    grad = partial(
        linear_loss.gradient,
        X=X[:, :-1],
        y=y,
        sample_weight=sw,
        l2_reg_strength=l2_reg_strength,
    )
    coef_penalized_without_intercept = _special_minimize(
        fun, grad, coef_unpenalized[:-1], tol_NM=1e-6, tol=1e-14
    )

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
        l2_reg_strength,
    )


@pytest.mark.parametrize("solver", SOLVERS)
@pytest.mark.parametrize("fit_intercept", [False, True])
def test_glm_regression(solver, fit_intercept, glm_dataset):
    """Test that GLM converges for all solvers to correct solution.

    We work with a simple constructed data set with known solution.
    """
    model, X, y, _, coef_with_intercept, coef_without_intercept, alpha = glm_dataset
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

    if solver == "lbfgs":
        rtol = 5e-5
    elif solver == "newton-lsmr":
        rtol = 5e-9
    else:
        rtol = 1e-9
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
    For long X, [X, X] is still a long but singular matrix.
    """
    model, X, y, _, coef_with_intercept, coef_without_intercept, alpha = glm_dataset
    n_samples, n_features = X.shape
    params = dict(
        alpha=alpha / 2,
        fit_intercept=fit_intercept,
        solver=solver,
        tol=1e-12,
        max_iter=1000,
    )

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

    with warnings.catch_warnings():
        # XXX: Investigate if the ConvergenceWarning that can appear in some
        # cases should be considered a bug or not. In the mean time we don't
        # fail when the assertions below pass irrespective of the presence of
        # the warning.
        warnings.simplefilter("ignore", ConvergenceWarning)
        model.fit(X, y)

    if solver == "lbfgs":
        rtol = 2e-4
    elif solver == "newton-lsmr":
        rtol = 2e-8
    else:
        rtol = 5e-9
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
    model, X, y, _, coef_with_intercept, coef_without_intercept, alpha = glm_dataset
    n_samples, n_features = X.shape
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

    rtol = 3e-5 if solver == "lbfgs" else 5e-9
    assert model.intercept_ == pytest.approx(intercept, rel=rtol)
    assert_allclose(model.coef_, coef, rtol=rtol)


@pytest.mark.parametrize("solver", SOLVERS)
@pytest.mark.parametrize("fit_intercept", [True, False])
def test_glm_regression_unpenalized(solver, fit_intercept, glm_dataset):
    """Test that unpenalized GLM converges for all solvers to correct solution.

    We work with a simple constructed data set with known solution.
    Note: This checks the minimum norm solution for wide X, i.e.
    n_samples < n_features:
        min ||w||_2 subject to w = argmin deviance(X, y, w)
    """
    model, X, y, coef, _, _, _ = glm_dataset
    n_samples, n_features = X.shape
    alpha = 0  # unpenalized
    params = dict(
        alpha=alpha,
        fit_intercept=fit_intercept,
        solver=solver,
        tol=1e-12,
        max_iter=1000,
    )

    model = clone(model).set_params(**params)
    if fit_intercept:
        X = X[:, :-1]  # remove intercept
        intercept = coef[-1]
        coef = coef[:-1]
    else:
        intercept = 0

    with warnings.catch_warnings():
        if solver.startswith("newton") and n_samples < n_features:
            # The newton solvers should warn and automatically fallback to LBFGS
            # in this case. The model should still converge.
            warnings.filterwarnings("ignore", category=scipy.linalg.LinAlgWarning)
        # XXX: Investigate if the ConvergenceWarning that can appear in some
        # cases should be considered a bug or not. In the mean time we don't
        # fail when the assertions below pass irrespective of the presence of
        # the warning.
        warnings.filterwarnings("ignore", category=ConvergenceWarning)
        model.fit(X, y)

    # FIXME: `assert_allclose(model.coef_, coef)` should work for all cases but fails
    # for the wide/fat case with n_features > n_samples. Most current GLM solvers do
    # NOT return the minimum norm solution with fit_intercept=True.
    if n_samples > n_features:
        rtol = 5e-5 if solver == "lbfgs" else 1e-7
        assert model.intercept_ == pytest.approx(intercept)
        assert_allclose(model.coef_, coef, rtol=rtol)
    else:
        # As it is an underdetermined problem, prediction = y. The following shows that
        # we get a solution, i.e. a (non-unique) minimum of the objective function ...
        rtol = 5e-5
        if solver == "newton-cholesky":
            rtol = 5e-4
        assert_allclose(model.predict(X), y, rtol=rtol)

        norm_solution = np.linalg.norm(np.r_[intercept, coef])
        norm_model = np.linalg.norm(np.r_[model.intercept_, model.coef_])
        if solver in ("newton-cholesky", "newton-lsmr"):
            # XXX: This solver shows random behaviour. Sometimes it finds solutions
            # with norm_model <= norm_solution! So we check conditionally.
            if norm_model < (1 + 1e-12) * norm_solution:
                assert model.intercept_ == pytest.approx(intercept)
                assert_allclose(model.coef_, coef, rtol=rtol)
        elif solver == "lbfgs" and fit_intercept:
            # But it is not the minimum norm solution. Otherwise the norms would be
            # equal.
            assert norm_model > (1 + 1e-12) * norm_solution

            # See https://github.com/scikit-learn/scikit-learn/issues/23670.
            # Note: Even adding a tiny penalty does not give the minimal norm solution.
            # XXX: We could have naively expected LBFGS to find the minimal norm
            # solution by adding a very small penalty. Even that fails for a reason we
            # do not properly understand at this point.
        else:
            # When `fit_intercept=False`, LBFGS naturally converges to the minimum norm
            # solution on this problem.
            # XXX: Do we have any theoretical guarantees why this should be the case?
            assert model.intercept_ == pytest.approx(intercept, rel=rtol)
            assert_allclose(model.coef_, coef, rtol=rtol)


@pytest.mark.parametrize("solver", SOLVERS)
@pytest.mark.parametrize("fit_intercept", [True, False])
def test_glm_regression_unpenalized_hstacked_X(solver, fit_intercept, glm_dataset):
    """Test that unpenalized GLM converges for all solvers to correct solution.

    We work with a simple constructed data set with known solution.
    GLM fit on [X] is the same as fit on [X, X]/2.
    For long X, [X, X] is a singular matrix and we check against the minimum norm
    solution:
        min ||w||_2 subject to w = argmin deviance(X, y, w)
    """
    model, X, y, coef, _, _, _ = glm_dataset
    n_samples, n_features = X.shape
    alpha = 0  # unpenalized
    params = dict(
        alpha=alpha,
        fit_intercept=fit_intercept,
        solver=solver,
        tol=1e-12,
        max_iter=1000,
    )

    model = clone(model).set_params(**params)
    if fit_intercept:
        intercept = coef[-1]
        coef = coef[:-1]
        if n_samples > n_features:
            X = X[:, :-1]  # remove intercept
            X = 0.5 * np.concatenate((X, X), axis=1)
        else:
            # To know the minimum norm solution, we keep one intercept column and do
            # not divide by 2. Later on, we must take special care.
            X = np.c_[X[:, :-1], X[:, :-1], X[:, -1]]
    else:
        intercept = 0
        X = 0.5 * np.concatenate((X, X), axis=1)
    assert np.linalg.matrix_rank(X) <= min(n_samples, n_features)

    with warnings.catch_warnings():
        if solver.startswith("newton"):
            # The newton solvers should warn and automatically fallback to LBFGS
            # in this case. The model should still converge.
            warnings.filterwarnings("ignore", category=scipy.linalg.LinAlgWarning)
        # XXX: Investigate if the ConvergenceWarning that can appear in some
        # cases should be considered a bug or not. In the mean time we don't
        # fail when the assertions below pass irrespective of the presence of
        # the warning.
        warnings.filterwarnings("ignore", category=ConvergenceWarning)
        model.fit(X, y)

    if fit_intercept and n_samples < n_features:
        # Here we take special care.
        model_intercept = 2 * model.intercept_
        model_coef = 2 * model.coef_[:-1]  # exclude the other intercept term.
        # For minimum norm solution, we would have
        # assert model.intercept_ == pytest.approx(model.coef_[-1])
    else:
        model_intercept = model.intercept_
        model_coef = model.coef_

    if n_samples > n_features:
        assert model_intercept == pytest.approx(intercept)
        rtol = 1e-4
        assert_allclose(model_coef, np.r_[coef, coef], rtol=rtol)
    else:
        # As it is an underdetermined problem, prediction = y. The following shows that
        # we get a solution, i.e. a (non-unique) minimum of the objective function ...
        rtol = 1e-6 if solver == "lbfgs" else 5e-6
        assert_allclose(model.predict(X), y, rtol=rtol)
        if (solver == "lbfgs" and fit_intercept) or solver in (
            "newton-cholesky",
            "newton-lsmr",
        ):
            # Same as in test_glm_regression_unpenalized.
            # But it is not the minimum norm solution. Otherwise the norms would be
            # equal.
            norm_solution = np.linalg.norm(
                0.5 * np.r_[intercept, intercept, coef, coef]
            )
            norm_model = np.linalg.norm(np.r_[model.intercept_, model.coef_])
            assert norm_model > (1 + 1e-12) * norm_solution
            # For minimum norm solution, we would have
            # assert model.intercept_ == pytest.approx(model.coef_[-1])
        else:
            assert model_intercept == pytest.approx(intercept, rel=5e-6)
            assert_allclose(model_coef, np.r_[coef, coef], rtol=1e-4)


@pytest.mark.parametrize("solver", SOLVERS)
@pytest.mark.parametrize("fit_intercept", [True, False])
def test_glm_regression_unpenalized_vstacked_X(solver, fit_intercept, glm_dataset):
    """Test that unpenalized GLM converges for all solvers to correct solution.

    We work with a simple constructed data set with known solution.
    GLM fit on [X] is the same as fit on [X], [y]
                                         [X], [y].
    For wide X, [X', X'] is a singular matrix and we check against the minimum norm
    solution:
        min ||w||_2 subject to w = argmin deviance(X, y, w)
    """
    model, X, y, coef, _, _, _ = glm_dataset
    n_samples, n_features = X.shape
    alpha = 0  # unpenalized
    params = dict(
        alpha=alpha,
        fit_intercept=fit_intercept,
        solver=solver,
        tol=1e-12,
        max_iter=1000,
    )

    model = clone(model).set_params(**params)
    if fit_intercept:
        X = X[:, :-1]  # remove intercept
        intercept = coef[-1]
        coef = coef[:-1]
    else:
        intercept = 0
    X = np.concatenate((X, X), axis=0)
    assert np.linalg.matrix_rank(X) <= min(n_samples, n_features)
    y = np.r_[y, y]

    with warnings.catch_warnings():
        if solver.startswith("newton") and n_samples < n_features:
            # The newton solvers should warn and automatically fallback to LBFGS
            # in this case. The model should still converge.
            warnings.filterwarnings("ignore", category=scipy.linalg.LinAlgWarning)
        # XXX: Investigate if the ConvergenceWarning that can appear in some
        # cases should be considered a bug or not. In the mean time we don't
        # fail when the assertions below pass irrespective of the presence of
        # the warning.
        warnings.filterwarnings("ignore", category=ConvergenceWarning)
        model.fit(X, y)

    if n_samples > n_features:
        rtol = 5e-5 if solver == "lbfgs" else 1e-6
        assert model.intercept_ == pytest.approx(intercept)
        assert_allclose(model.coef_, coef, rtol=rtol)
    else:
        # As it is an underdetermined problem, prediction = y. The following shows that
        # we get a solution, i.e. a (non-unique) minimum of the objective function ...
        rtol = 1e-6 if solver == "lbfgs" else 5e-6
        assert_allclose(model.predict(X), y, rtol=rtol)

        norm_solution = np.linalg.norm(np.r_[intercept, coef])
        norm_model = np.linalg.norm(np.r_[model.intercept_, model.coef_])
        if solver in ("newton-cholesky", "newton-lsmr"):
            # XXX: This solver shows random behaviour. Sometimes it finds solutions
            # with norm_model <= norm_solution! So we check conditionally.
            if not (norm_model > (1 + 1e-12) * norm_solution):
                assert model.intercept_ == pytest.approx(intercept)
                assert_allclose(model.coef_, coef, rtol=1e-4)
        elif solver == "lbfgs" and fit_intercept:
            # Same as in test_glm_regression_unpenalized.
            # But it is not the minimum norm solution. Otherwise the norms would be
            # equal.
            assert norm_model > (1 + 1e-12) * norm_solution
        else:
            rtol = 1e-5 if solver in ("newton-cholesky", "newton-lsmr") else 1e-4
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
    alpha = 1
    params = {
        "solver": solver,
        "fit_intercept": fit_intercept,
        "tol": 1e-10,
    }

    glm1 = PoissonRegressor(warm_start=False, max_iter=1000, alpha=alpha, **params)
    glm1.fit(X, y)

    glm2 = PoissonRegressor(warm_start=True, max_iter=1, alpha=alpha, **params)
    # As we intentionally set max_iter=1 such that the solver should raise a
    # ConvergenceWarning.
    with pytest.warns(ConvergenceWarning):
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
        l2_reg_strength=alpha,
    )
    objective_glm2 = linear_loss.loss(
        coef=np.r_[glm2.coef_, glm2.intercept_] if fit_intercept else glm2.coef_,
        X=X,
        y=y,
        sample_weight=sw,
        l2_reg_strength=alpha,
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


@pytest.mark.parametrize("solver", ["lbfgs", "newton-cholesky", "newton-lsmr"])
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


@pytest.mark.parametrize("newton_solver", ["newton-cholesky", "newton-lsmr"])
def test_linalg_warning_with_newton_solver(newton_solver, global_random_seed):
    rng = np.random.RandomState(global_random_seed)
    # Use at least 20 samples to reduce the likelihood of getting a degenerate
    # dataset for any global_random_seed.
    X_orig = rng.normal(size=(20, 3))
    y = rng.poisson(
        np.exp(X_orig @ np.ones(X_orig.shape[1])), size=X_orig.shape[0]
    ).astype(np.float64)

    # Collinear variation of the same input features.
    X_collinear = np.hstack([X_orig] * 10)

    # Let's consider the deviance of a constant baseline on this problem.
    baseline_pred = np.full_like(y, y.mean())
    constant_model_deviance = mean_poisson_deviance(y, baseline_pred)
    assert constant_model_deviance > 1.0

    # No warning raised on well-conditioned design, even without regularization.
    tol = 1e-10
    with warnings.catch_warnings():
        warnings.simplefilter("error")
        reg = PoissonRegressor(solver=newton_solver, alpha=0.0, tol=tol).fit(X_orig, y)
    original_newton_deviance = mean_poisson_deviance(y, reg.predict(X_orig))

    # On this dataset, we should have enough data points to not make it
    # possible to get a near zero deviance (for the any of the admissible
    # random seeds). This will make it easier to interpret meaning of rtol in
    # the subsequent assertions:
    assert original_newton_deviance > 0.2

    # We check that the model could successfully fit information in X_orig to
    # improve upon the constant baseline by a large margin (when evaluated on
    # the traing set).
    assert constant_model_deviance - original_newton_deviance > 0.1

    # LBFGS is robust to a collinear design because its approximation of the
    # Hessian is Symmeric Positive Definite by construction. Let's record its
    # solution
    with warnings.catch_warnings():
        warnings.simplefilter("error")
        reg = PoissonRegressor(solver="lbfgs", alpha=0.0, tol=tol).fit(X_collinear, y)
    collinear_lbfgs_deviance = mean_poisson_deviance(y, reg.predict(X_collinear))

    # The LBFGS solution on the collinear is expected to reach a comparable
    # solution to the Newton solution on the original data.
    rtol = 1e-6
    assert collinear_lbfgs_deviance == pytest.approx(original_newton_deviance, rel=rtol)

    # Fitting a Newton solver on the collinear version of the training data
    # without regularization should raise an informative warning and fallback
    # to the LBFGS solver.
    if newton_solver == "newton-cholesky":
        msg = (
            "The inner solver of .*Newton.*Solver stumbled upon a singular or very "
            "ill-conditioned Hessian matrix"
        )
        with pytest.warns(scipy.linalg.LinAlgWarning, match=msg):
            reg = PoissonRegressor(solver=newton_solver, alpha=0.0, tol=tol).fit(
                X_collinear, y
            )
    else:
        # newton-lsmr has no problems with collinearity
        reg = PoissonRegressor(solver=newton_solver, alpha=0.0, tol=tol).fit(
            X_collinear, y
        )
    # As a result we should still automatically converge to a good solution.
    collinear_newton_deviance = mean_poisson_deviance(y, reg.predict(X_collinear))
    assert collinear_newton_deviance == pytest.approx(
        original_newton_deviance, rel=rtol
    )

    # Increasing the regularization slightly should make the problem go away:
    with warnings.catch_warnings():
        warnings.simplefilter("error", scipy.linalg.LinAlgWarning)
        reg = PoissonRegressor(solver=newton_solver, alpha=1e-10).fit(X_collinear, y)

    # The slightly penalized model on the collinear data should be close enough
    # to the unpenalized model on the original data.
    penalized_collinear_newton_deviance = mean_poisson_deviance(
        y, reg.predict(X_collinear)
    )
    assert penalized_collinear_newton_deviance == pytest.approx(
        original_newton_deviance, rel=rtol
    )


@pytest.mark.parametrize(
    "solver, warn1, msg1, warn2, msg2",
    [
        (
            "lbfgs",
            None,
            None,
            [RuntimeWarning],
            ["invalid value encountered in matmul"],
        ),
        (
            "newton-cholesky",
            scipy.linalg.LinAlgWarning,
            (
                "The inner solver of .*Newton.*Solver stumbled upon a singular or very "
                "ill-conditioned Hessian matrix"
            ),
            [scipy.linalg.LinAlgWarning, RuntimeWarning],
            [
                (
                    "The inner solver of .*Newton.*Solver stumbled upon a singular or"
                    " very ill-conditioned Hessian matrix"
                ),
                "invalid value encountered in matmul",
            ],
        ),
        (
            "newton-lsmr",
            None,
            None,
            [ConvergenceWarning],
            ["Newton solver did not converge after .* iterations"],
        ),
    ],
)
def test_solver_on_ill_conditioned_X(
    solver, warn1, msg1, warn2, msg2, global_random_seed
):
    """Test GLM solvers on ill conditioned X with high condition number.

    Note that numbers in this test are tuned such that is passes for all global random
    seeds.
    """
    rng = np.random.RandomState(global_random_seed)
    # Use at least 20 samples to reduce the likelihood of getting a degenerate
    # dataset for any global_random_seed.
    X_orig = rng.uniform(low=-1, high=1, size=(20, 3))
    y = rng.poisson(
        np.exp(X_orig @ np.ones(X_orig.shape[1])), size=X_orig.shape[0]
    ).astype(np.float64)

    tol = 1e-7
    model = PoissonRegressor(solver=solver, alpha=0.0, tol=tol)

    # No warning raised on well-conditioned design, even without regularization.
    with warnings.catch_warnings():
        warnings.simplefilter("error")
        reg = clone(model).fit(X_orig, y)
    original_deviance = mean_poisson_deviance(y, reg.predict(X_orig))

    # Construct an ill-conditioned problem by some almost identical columns.
    X_ill_conditioned = X_orig.copy()
    X_ill_conditioned[:, 1] = X_ill_conditioned[:, 0]
    X_ill_conditioned[0, 1] += 1e-10
    X_ill_conditioned[-1, 1] -= 1e-10
    # Make sure that it is ill-conditioned <=> large condition number.
    assert np.linalg.cond(X_ill_conditioned) > 1e10 * np.linalg.cond(X_orig)

    if warn1 is None:
        # Assert that no warning is raised.
        with warnings.catch_warnings():
            warnings.simplefilter("error")
            reg = clone(model).fit(X_ill_conditioned, y)
    else:
        # Assert that the given warning is raised.
        with pytest.warns(warn1, match=msg1):
            reg = clone(model).fit(X_ill_conditioned, y)

    # Construct another ill-conditioned problem by scaling of features.
    X_ill_conditioned = X_orig.copy()
    X_ill_conditioned[:, 0] *= 1e-4
    X_ill_conditioned[:, 1] *= 1e4
    # Make sure that it is ill conditioned >=> large condition number.
    assert np.linalg.cond(X_ill_conditioned) > 1e5 * np.linalg.cond(X_orig)

    test_loss = False
    if warn2 is None:
        with warnings.catch_warnings():
            warnings.simplefilter("error")
            reg = clone(model).fit(X_ill_conditioned, y)
    else:
        # Whether or not a warning is raised depends on global_random_seed.
        with warnings.catch_warnings(record=True) as w:
            for warn in warn2:
                warnings.simplefilter("ignore", warn)
            reg = clone(model).fit(X_ill_conditioned, y)
            test_loss = True
            assert len(w) == 0

    if test_loss:
        # Without penalty, scaling of columns has no effect on predictions.
        ill_cond_deviance = mean_poisson_deviance(y, reg.predict(X_ill_conditioned))
        if solver in ("lbfgs", "newton-cholesky", "newton-lsmr"):
            pytest.xfail(
                f"Solver {solver} does not converge but does so without warning."
            )
        assert ill_cond_deviance == pytest.approx(original_deviance, rel=1e-2)


@pytest.mark.parametrize("verbose", [0, 1, 2])
def test_newton_solver_verbosity(capsys, verbose):
    """Test the std output of verbose newton solvers."""
    y = np.array([1, 2], dtype=float)
    X = np.array([[1.0, 0], [0, 1]], dtype=float)
    linear_loss = LinearModelLoss(base_loss=HalfPoissonLoss(), fit_intercept=False)
    sol = NewtonCholeskySolver(
        coef=linear_loss.init_zero_coef(X),
        linear_loss=linear_loss,
        l2_reg_strength=0,
        verbose=verbose,
    )
    sol.solve(X, y, None)  # returns array([0., 0.69314758])
    captured = capsys.readouterr()

    if verbose == 0:
        assert captured.out == ""
    else:
        msg = [
            "Newton iter=1",
            "Check Convergence",
            "1. max |gradient|",
            "2. Newton decrement",
            "Solver did converge at loss = ",
        ]
        for m in msg:
            assert m in captured.out

    if verbose >= 2:
        msg = ["Backtracking Line Search", "line search iteration="]
        for m in msg:
            assert m in captured.out

    # Set the Newton solver to a state with a completely wrong Newton step.
    sol = NewtonCholeskySolver(
        coef=linear_loss.init_zero_coef(X),
        linear_loss=linear_loss,
        l2_reg_strength=0,
        verbose=verbose,
    )
    sol.setup(X=X, y=y, sample_weight=None)
    sol.iteration = 1
    sol.update_gradient_hessian(X=X, y=y, sample_weight=None)
    sol.coef_newton = np.array([1.0, 0])
    sol.gradient_times_newton = sol.gradient @ sol.coef_newton
    with warnings.catch_warnings():
        warnings.simplefilter("ignore", ConvergenceWarning)
        sol.line_search(X=X, y=y, sample_weight=None)
        captured = capsys.readouterr()
    if verbose >= 1:
        assert (
            "Line search did not converge and resorts to lbfgs instead." in captured.out
        )

    # Set the Newton solver to a state with bad Newton step such that the loss
    # improvement in line search is tiny.
    sol = NewtonCholeskySolver(
        coef=np.array([1e-12, 0.69314758]),
        linear_loss=linear_loss,
        l2_reg_strength=0,
        verbose=verbose,
    )
    sol.setup(X=X, y=y, sample_weight=None)
    sol.iteration = 1
    sol.update_gradient_hessian(X=X, y=y, sample_weight=None)
    sol.coef_newton = np.array([1e-6, 0])
    sol.gradient_times_newton = sol.gradient @ sol.coef_newton
    with warnings.catch_warnings():
        warnings.simplefilter("ignore", ConvergenceWarning)
        sol.line_search(X=X, y=y, sample_weight=None)
        captured = capsys.readouterr()
    if verbose >= 2:
        msg = [
            "line search iteration=",
            "check loss improvement <= armijo term:",
            "check loss |improvement| <= eps * |loss_old|:",
            "check sum(|gradient|) < sum(|gradient_old|):",
        ]
        for m in msg:
            assert m in captured.out

    # Test for a case with negative hessian. We badly initialize coef for a Tweedie
    # loss with non-canonical link, e.g. Inverse Gaussian deviance with a log link.
    linear_loss = LinearModelLoss(
        base_loss=HalfTweedieLoss(power=3), fit_intercept=False
    )
    sol = NewtonCholeskySolver(
        coef=linear_loss.init_zero_coef(X) + 1,
        linear_loss=linear_loss,
        l2_reg_strength=0,
        verbose=verbose,
    )
    with warnings.catch_warnings():
        warnings.simplefilter("ignore", ConvergenceWarning)
        sol.solve(X, y, None)
    captured = capsys.readouterr()
    if verbose >= 1:
        assert (
            "The inner solver detected a pointwise Hessian with many negative values"
            " and resorts to lbfgs instead."
            in captured.out
        )


@pytest.mark.skipif(
    sp_version < parse_version("1.4.0"),
    reason="LinearOperator transpose needs scipy>=1.4.0",
)
@pytest.mark.parametrize("fit_intercept", [False, True])
@pytest.mark.parametrize("l2_reg_strength", [0, 1.5])
def test_NewtonLSMRSolver_multinomial_A_b(
    fit_intercept, l2_reg_strength, global_random_seed
):
    """Test NewtonLSMRSolver.compute_A_b for multinomial case"""
    n_samples, n_features, n_classes = 12, 10, 5
    n_dof = n_features + fit_intercept
    rng = np.random.RandomState(global_random_seed)
    X, y = make_classification(
        n_samples=n_samples,
        n_features=n_features,
        n_classes=n_classes,
        n_informative=n_features - 1,
        n_redundant=1,
        random_state=rng,
    )
    y = y.astype(float)
    coef = rng.standard_normal(size=(n_classes, n_features + fit_intercept))
    coef -= np.mean(coef, axis=0)

    multinomial_loss = LinearModelLoss(
        base_loss=HalfMultinomialLoss(n_classes=n_classes), fit_intercept=fit_intercept
    )
    gradient, hessp = multinomial_loss.gradient_hessian_product(
        coef=coef,
        X=X,
        y=y,
        sample_weight=None,
        l2_reg_strength=l2_reg_strength,
    )
    sol = NewtonLSMRSolver(
        coef=coef,
        linear_loss=multinomial_loss,
        l2_reg_strength=l2_reg_strength,
    )
    sol.setup(X=X, y=y, sample_weight=None)
    sol.update_gradient_hessian(X=X, y=y, sample_weight=None)

    # Hessian @ coef
    # with hessp
    H_coef = hessp(coef)  # shape = (n_classes, n_dof)
    # with A
    A, b = sol.compute_A_b(X=X, y=y, sample_weight=None)
    At_A_coef = A.T @ A @ coef.ravel(order="C")
    # The results are differently ravelled, we restore the 2-d arrays with n_classes
    # on the 1st (=last) axis.
    assert_allclose(At_A_coef.reshape(-1, n_classes, order="F"), H_coef.T)
    # Note: The following does not work for all global_random_seeds. The reason behind
    # it is unclear.
    # assert_allclose(
    #     A.rmatvec(b).reshape(-1, n_classes, order="F"), -gradient.T, rtol=1e-4
    # )

    # Test consistency of A, i.e. reconstructing the matrix based on
    # A @ unit_vector and A.T @ unit_vector should give the same matrix.
    unit_vec = np.zeros(n_dof * n_classes)
    A_matrix1 = np.zeros(((n_samples + n_features) * n_classes, n_dof * n_classes))
    for i in range(n_dof * n_classes):
        unit_vec[i] = 1
        A_matrix1[:, i] = A @ unit_vec
        unit_vec[i] = 0
    unit_vec = np.zeros((n_samples + n_features) * n_classes)
    A_matrix2 = np.zeros(((n_samples + n_features) * n_classes, n_dof * n_classes))
    for j in range((n_samples + n_features) * n_classes):
        unit_vec[j] = 1
        A_matrix2[j, :] = A.rmatvec(unit_vec)
        unit_vec[j] = 0
    assert_allclose(A_matrix1, A_matrix2)


@pytest.mark.skipif(
    sp_version < parse_version("1.4.0"),
    reason="LinearOperator transpose needs scipy>=1.4.0",
)
@pytest.mark.parametrize("fit_intercept", (False, True))
@pytest.mark.parametrize("with_sample_weight", (False, True))
def test_NewtonLSMRSolver_multinomial_A_b_on_3_classes(
    fit_intercept, with_sample_weight, global_random_seed
):
    """Test NewtonLSMRSolver.compute_A_b for multinomial case with 3 classes."""
    n_samples, n_features, n_classes = 5, 2, 3
    n_dof = n_features + fit_intercept
    rng = np.random.RandomState(global_random_seed)
    X = np.array([[1.0, 0], [1, 1], [1, 2], [1, 3], [1, 4]])
    Xi = X
    # coef.shape = (n_classes, n_dof) and coef.sum(axis=0) = 0
    coef = np.array([[1.0, -1], [1, 0], [-2, 1]])
    if fit_intercept:
        X[0, 0] = -1  # to make Xi positive definite
        coef = np.c_[coef, [[2], [1], [-3]]]
        Xi = np.c_[X, np.ones(n_samples)]
    assert np.linalg.matrix_rank(Xi) == n_dof
    y = np.array([0.0, 1, 1, 1, 2])

    if with_sample_weight:
        sw = 1.0 + np.arange(n_samples)
    else:
        sw = None

    multinomial_loss = LinearModelLoss(
        base_loss=HalfMultinomialLoss(n_classes=n_classes, sample_weight=sw),
        fit_intercept=fit_intercept,
    )
    sol = NewtonLSMRSolver(
        coef=coef,
        linear_loss=multinomial_loss,
        l2_reg_strength=0,
    )
    sol.setup(X=X, y=y, sample_weight=sw)
    sol.update_gradient_hessian(X=X, y=y, sample_weight=sw)
    A, b = sol.compute_A_b(X=X, y=y, sample_weight=sw)

    _, _, raw_prediction = multinomial_loss.weight_intercept_raw(coef, X)
    # p.shape = (n_samples, n_classes)
    p = multinomial_loss.base_loss.predict_proba(raw_prediction)
    assert_allclose(np.sum(p, axis=1), 1)  # consequence of coef.sum(axis=0) = 0

    # Extended vectors for 3 classes:
    # For y_ext and p_ext, we append the vectors for different classes one after the
    # other.
    y_ext = np.r_[y == 0, y == 1, y == 2]
    p_ext = p.ravel(order="F")
    coef_ext = coef.ravel(order="C")
    sw_ext = np.tile(sw, 3)
    # Extended matrices for 3 classes:
    # X_ext = [X, 0, 0]  W = [W00, W01, W02]
    #         [0, X, 0]      [W10, W11, W12]
    #         [0, 0, X]      [W20, W21, W22]
    # X_ext.shape = (n_classes * n_samples, n_classes * n_features)
    # W.shape = (n_classes * n_samples, n_classes * n_samples)
    # Each Wij is a diagonal matrix with elements pi * (delta_ij - pj).
    X_ext = linalg.block_diag(Xi, Xi, Xi)
    W00 = np.diag(p[:, 0] * (1 - p[:, 0]))
    W11 = np.diag(p[:, 1] * (1 - p[:, 1]))
    W22 = np.diag(p[:, 2] * (1 - p[:, 2]))
    W01 = np.diag(-p[:, 0] * p[:, 1])
    W02 = np.diag(-p[:, 0] * p[:, 2])
    W12 = np.diag(-p[:, 1] * p[:, 2])
    W_ext = np.block([[W00, W01, W02], [W01, W11, W12], [W02, W12, W22]])
    g_ext = p_ext - y_ext  # pointwise gradient
    if with_sample_weight:
        W_ext *= sw_ext[:, None]
        g_ext *= sw_ext
    G_ext = X_ext.T @ g_ext  # full gradient
    H_ext = X_ext.T @ W_ext @ X_ext  # full hessian

    assert_allclose(g_ext.reshape((n_samples, n_classes), order="F"), sol.g)

    # Note that both X_ext @ v and A @ v need v = coef.ravel(order="C") for
    # coef.shape = (n_classes, n_dof).

    # v is like a coefficient
    v = rng.standard_normal(size=(n_classes, n_dof))
    # u is a test vector for A', i.e. we test A.T @ u
    u = rng.standard_normal(size=(n_samples + n_features, n_classes))

    # Check that A is a square root of H, i.e. H = A'A.
    H_ext_v = H_ext @ v.ravel(order="C")
    A_v = A @ v.ravel(order="C")
    At_A_v = A.T @ A_v
    assert_allclose(At_A_v, H_ext_v)

    # Check that A'b = -G
    assert_allclose(A.T @ b, -G_ext)

    # Check that A is constructed from LDL decomposition.
    LDL = Multinomial_LDL_Decomposition(proba=p)
    D_Lt_X_v = LDL.sqrt_D_Lt_matmul(Xi @ v.T)
    if with_sample_weight:
        D_Lt_X_v *= np.sqrt(sw)[:, None]
    assert_allclose(A_v.reshape((-1, n_classes), order="F")[:n_samples, :], D_Lt_X_v)
    assert_allclose(A_v[n_classes * n_samples :], 0)
    assert_allclose(
        A_v.reshape((-1, n_classes), order="F")[:n_samples, :], D_Lt_X_v, atol=1e-8
    )

    # Check LDL decomposition of W explicitly.
    L, D, perm = linalg.ldl(W_ext, lower=True)
    # D should be diagonal (and not block diagonal) and non-negative.
    assert_allclose(np.diag(np.diag(D)), D)
    assert np.all(D >= -1e-15)
    D[D < 1e-15] = 0  # Tiny values would be zero in exact arithmetic.
    assert_allclose(L @ D @ L.T, W_ext)
    D_Lt_X = np.sqrt(D) @ L.T @ X_ext
    assert_allclose(D_Lt_X.T @ D_Lt_X, H_ext)
    assert_allclose(D_Lt_X @ v.ravel(order="C"), D_Lt_X_v.ravel(order="F"), atol=1e-8)

    # Construct matrix A explicitly. Note: no intercept, penalty is zero.
    A_ext = np.r_[D_Lt_X, np.zeros((n_classes * n_dof, n_classes * n_dof))]
    assert_allclose(A_ext.T @ A_ext, H_ext)
    D_inv = np.sqrt(D)
    D_inv[D_inv > 0] = 1 / D_inv[D_inv > 0]
    D_inv[D_inv <= 0] = 0
    assert_allclose(L @ np.sqrt(D) @ D_inv @ linalg.solve(L, g_ext), g_ext)
    assert_allclose(X_ext.T @ L @ np.sqrt(D) @ D_inv @ linalg.solve(L, g_ext), G_ext)
    # Note that A_ext needs slicing u in a way not possible with ravel, i.e.
    # treat u[:n_samples * n_classes, :] and u[n_samples * n_classes:, :] differently.
    assert_allclose(
        A_ext[: n_classes * n_samples, :].T @ u[:n_samples, :].ravel(order="F"),
        A.T @ u.ravel(order="F"),
    )

    # Construct b explicitly.
    b_ext = np.r_[-D_inv @ linalg.solve(L, g_ext), np.zeros(n_dof * n_classes)]
    assert_allclose(
        b_ext[: n_samples * n_classes],
        b.reshape((-1, n_classes), order="F")[:n_samples, :].ravel(order="F"),
        atol=1e-15,
    )
    assert_allclose(A_ext.T @ b_ext, -G_ext)

    # Check equivalence of linear equation H @ x = H_v to least squares problem
    # ||D_Lt_X @ x - D_Lt_X_v||. Note that H (unpenalized with all classes) is
    # singular and we use `lstsq` instead of `solve`.
    res = linalg.lstsq(H_ext, H_ext_v)[0].reshape((n_classes, -1))
    assert_allclose(res, v - np.mean(v, axis=0))
    res = linalg.lstsq(D_Lt_X, D_Lt_X_v.ravel(order="F"))[0].reshape((n_classes, -1))
    assert_allclose(res - np.mean(res, axis=0), v - np.mean(v, axis=0))

    # Check Newton step, i.e. H @ x = -G.
    # Here, we need to add a penalty term. Otherwise, the linear system would be
    # singular.
    alpha = 1e-4
    sol = NewtonLSMRSolver(
        coef=coef,
        linear_loss=multinomial_loss,
        l2_reg_strength=alpha,
    )
    sol.setup(X=X, y=y, sample_weight=sw)
    sol.update_gradient_hessian(X=X, y=y, sample_weight=sw)
    A, b = sol.compute_A_b(X=X, y=y, sample_weight=sw)
    pen_ext = alpha * np.identity(n_classes * n_dof)
    if fit_intercept:
        pen_ext[:, n_dof - 1 :: n_dof] = 0
    assert_allclose(
        A.T @ A @ v.ravel(order="C"), (H_ext + pen_ext) @ v.ravel(order="C")
    )
    assert_allclose(A.T @ b, -(G_ext + pen_ext @ coef_ext))

    # Check Newton step.
    # Note: This works despite the fact that W=LDL is a singular matrix and therefore
    # D has some zeros on the diagonal and is not strictly invertible. But this does
    # not seem to matter.
    if fit_intercept:
        with warnings.catch_warnings():
            # Due to the intercept term not being penalized, the linear system might
            # still be singular.
            warnings.simplefilter("ignore", scipy.linalg.LinAlgWarning)
            res1 = linalg.solve(H_ext + pen_ext, -(G_ext + pen_ext @ coef_ext)).reshape(
                (n_classes, -1)
            )
        # For the same reasons, the intercept might need class centering.
        res1[:, -1] -= np.mean(res1[:, -1])
    else:
        res1 = linalg.solve(H_ext + pen_ext, -(G_ext + pen_ext @ coef_ext)).reshape(
            (n_classes, -1)
        )
    # sum over classes = 0
    assert_allclose(res1.sum(axis=0), 0, atol=5e-9)
    assert_allclose((H_ext + pen_ext) @ res1.ravel(), -(G_ext + pen_ext @ coef_ext))
    assert_allclose(A.T @ A @ res1.ravel(order="C"), A.T @ b)
    res2 = lsmr(A, b, maxiter=(n_features + n_samples) * n_classes, show=True)
    res2 = res2[0].reshape((n_classes, -1))
    assert_allclose(res1, res2)


@pytest.mark.parametrize("fit_intercept", [False, True])
@pytest.mark.parametrize("l2_reg_strength", [0, 5.5])
def test_NewtonLSMRSolver_multinomial_on_binary_problem(
    fit_intercept, l2_reg_strength, global_random_seed
):
    """Test NewtonLSMRSolver multinomial case on a binary problem."""
    n_samples, n_features, n_classes = 100, 3, 2
    rng = np.random.RandomState(global_random_seed)
    X, y = make_classification(
        n_samples=n_samples,
        n_features=n_features,
        n_classes=n_classes,
        n_informative=n_features - 1 * 0,
        n_redundant=1 * 0,
        flip_y=0.1,
        random_state=rng,
    )
    y = y.astype(float)

    # Note that the factor of 2 for BinomialRegressor is to adjust for the
    # parametrization difference to the multinomial loss.
    # The factor 1/n_samples is due to the scaling happening in BinomialRegressor
    # in contrast to NewtonLSMRSolver.
    bin = BinomialRegressor(
        fit_intercept=fit_intercept,
        alpha=l2_reg_strength / 2 / n_samples,
        solver="newton-cholesky",
        tol=1e-8,
    ).fit(X, y)
    bin_loss = LinearModelLoss(
        base_loss=HalfBinomialLoss(),
        fit_intercept=fit_intercept,
    )

    linear_loss = LinearModelLoss(
        base_loss=HalfMultinomialLoss(n_classes=n_classes),
        fit_intercept=fit_intercept,
    )
    sol = NewtonLSMRSolver(
        coef=linear_loss.init_zero_coef(X),
        linear_loss=linear_loss,
        l2_reg_strength=l2_reg_strength,
        tol=1e-8,
        verbose=0,
    )
    sol.solve(X, y, None)

    assert_allclose(np.mean(sol.coef, axis=0), 0, atol=1e-13)
    if fit_intercept:
        coef_bin = np.r_[bin.coef_, bin.intercept_]
    else:
        coef_bin = bin.coef_
    assert_allclose(
        linear_loss.loss(coef=sol.coef, X=X, y=y, l2_reg_strength=l2_reg_strength),
        bin_loss.loss(coef=coef_bin, X=X, y=y, l2_reg_strength=l2_reg_strength / 2),
    )
    assert_allclose(sol.coef[1, :], coef_bin / 2, rtol=5e-6)
