# Authors: The scikit-learn developers
# SPDX-License-Identifier: BSD-3-Clause

import numpy as np
import pytest
from numpy.testing import assert_allclose, assert_array_equal
from scipy import sparse
from scipy.stats import norm

from sklearn.base import clone
from sklearn.linear_model import Lasso, LassoCV
from sklearn.model_selection import RepeatedKFold


@pytest.fixture
def data():
    rng = np.random.RandomState(12)
    X = rng.normal(size=(53, 6)) + 3
    return X, 2 + X[:, 0] + rng.normal(size=len(X))


def estimator(**kwargs):
    params = dict(
        alphas=np.geomspace(1, 0.01, 7),
        cv=RepeatedKFold(n_splits=3, n_repeats=4, random_state=12),
        cv_tolerance=100,
        min_cv_repeats=4,
        tol=1e-10,
        max_iter=10000,
    )
    params.update(kwargs)
    return LassoCV(**params)


@pytest.mark.parametrize("as_sparse", [False, True])
@pytest.mark.parametrize("weighted", [False, True])
@pytest.mark.parametrize("fit_intercept", [False, True])
@pytest.mark.parametrize("positive", [False, True])
def test_fixed_budget_equivalence(data, as_sparse, weighted, fit_intercept, positive):
    X, y = data
    if as_sparse:
        X = sparse.csc_matrix(X)
    weights = np.linspace(0.5, 2, len(y)) if weighted else None
    adaptive = estimator(fit_intercept=fit_intercept, positive=positive).fit(
        X, y, sample_weight=weights
    )
    fixed = (
        clone(adaptive).set_params(cv_tolerance=None).fit(X, y, sample_weight=weights)
    )
    assert_allclose(adaptive.mse_path_, fixed.mse_path_, rtol=1e-10, atol=1e-12)
    assert_allclose(adaptive.coef_, fixed.coef_, rtol=1e-10, atol=1e-12)
    assert adaptive.alpha_ == fixed.alpha_
    assert adaptive.n_repeats_ == 4
    assert adaptive.mse_path_.shape == (7, 12)
    curves = adaptive.mse_path_.reshape(7, 4, 3).mean(axis=2).T
    relative = curves - curves[:, :1]
    assert_allclose(adaptive.cv_mcse_, relative.std(axis=0, ddof=1) / 2, atol=1e-14)
    assert_allclose(
        adaptive.cv_reproducibility_threshold_, 100 / (np.sqrt(2) * norm.isf(0.05 / 12))
    )


def test_independent_fold_fits(data):
    X, y = data
    model = estimator().fit(X, y)
    expected = []
    for train, test in model.cv.split(X):
        expected.append(
            [
                np.mean(
                    (
                        y[test]
                        - Lasso(alpha=alpha, tol=1e-12, max_iter=10000)
                        .fit(X[train], y[train])
                        .predict(X[test])
                    )
                    ** 2
                )
                for alpha in model.alphas_
            ]
        )
    assert_allclose(model.mse_path_, np.asarray(expected).T, rtol=1e-8)


def test_work_count_and_reset(data, monkeypatch):
    calls, refits = [], []
    path, fit = LassoCV.path, Lasso.fit

    def counted_path(*args, **kwargs):
        calls.append(kwargs["alphas"].copy())
        return path(*args, **kwargs)

    def counted_fit(self, *args, **kwargs):
        refits.append(self.alpha)
        return fit(self, *args, **kwargs)

    monkeypatch.setattr(LassoCV, "path", staticmethod(counted_path))
    monkeypatch.setattr(Lasso, "fit", counted_fit)
    model = estimator(min_cv_repeats=2).fit(*data)
    assert len(calls) == 6
    assert len(refits) == 1
    assert model.n_repeats_ == 2
    assert model.cv_reproducibility_status_ == "target_reached"
    for grid in calls:
        assert_array_equal(grid, model.alphas_)
    model.set_params(cv_tolerance=None).fit(*data)
    for name in (
        "n_repeats_",
        "cv_mcse_",
        "cv_reproducibility_threshold_",
        "cv_reproducibility_status_",
    ):
        assert not hasattr(model, name)


def test_exhaustion_and_parallel(data):
    model = estimator(cv_tolerance=1e-15, min_cv_repeats=2)
    with pytest.warns(UserWarning, match="exhausted"):
        model.fit(*data)
    assert model.cv_reproducibility_status_ == "budget_exhausted"
    assert model.n_repeats_ == 4
    with pytest.warns(UserWarning, match="exhausted"):
        parallel = clone(model).set_params(n_jobs=2).fit(*data)
    assert_array_equal(model.mse_path_, parallel.mse_path_)
    assert_allclose(model.predict(data[0]), parallel.predict(data[0]))


@pytest.mark.parametrize(
    "kwargs",
    [
        {"cv": 3},
        {"min_cv_repeats": 5},
        {"min_cv_repeats": 1},
        {"cv_tolerance": 0},
        {"cv_tolerance": np.inf},
        {"cv_tolerance": np.nan},
        {"cv_error_rate": 0},
        {"cv_error_rate": 1},
        {"cv_error_rate": np.nan},
        {"alphas": [0, 1]},
        {"alphas": [1, 1]},
        {"alphas": [1]},
        {"selection": "random"},
        {"cv_error_rate": 5e-324},
    ],
)
def test_invalid_configuration(data, kwargs):
    with pytest.raises(ValueError):
        estimator(**kwargs).fit(*data)


def test_zero_variance(data):
    X, _ = data
    model = estimator(min_cv_repeats=2).fit(X, np.zeros(len(X)))
    assert model.n_repeats_ == 2
    assert_array_equal(model.cv_mcse_, 0)


def test_generated_grid_and_random_selection(data):
    model = estimator(alphas=9, selection="random", random_state=123).fit(*data)
    again = clone(model).fit(*data)
    assert model.alphas_.shape == (9,)
    assert_array_equal(model.mse_path_, again.mse_path_)


def test_feature_names(data):
    pd = pytest.importorskip("pandas")
    X, y = data
    frame = pd.DataFrame(X, columns=list("abcdef"))
    model = estimator().fit(frame, y)
    assert_array_equal(model.feature_names_in_, frame.columns)
    assert model.predict(frame).shape == y.shape
    with pytest.raises(ValueError, match="feature names"):
        model.predict(frame[frame.columns[::-1]])


@pytest.mark.parametrize("scale, expected_repeats", [(1, 3), (4, 10)])
def test_first_crossing_after_burn_in(data, monkeypatch, scale, expected_repeats):
    from sklearn.linear_model import _coordinate_descent

    differences = scale * np.array([0, 2] + [1] * 10)
    losses = iter(np.repeat(differences, 3))
    calls = []

    def residuals(*args, **kwargs):
        difference = next(losses)
        calls.append(difference)
        return np.array([10, 10 + difference])

    monkeypatch.setattr(_coordinate_descent, "_path_residuals", residuals)
    model = estimator(
        alphas=[1, 0.1],
        cv=RepeatedKFold(n_splits=3, n_repeats=12, random_state=12),
        min_cv_repeats=2,
        cv_tolerance=0.6 * np.sqrt(2) * norm.isf(0.025),
    ).fit(*data)
    assert model.n_repeats_ == expected_repeats
    assert len(calls) == 3 * expected_repeats
    assert_allclose(
        model.cv_mcse_[1],
        np.std(differences[:expected_repeats], ddof=1) / np.sqrt(expected_repeats),
    )
    assert model.cv_mcse_[1] <= 0.6
    assert (
        np.std(differences[: expected_repeats - 1], ddof=1)
        / np.sqrt(expected_repeats - 1)
        > 0.6
    )


@pytest.mark.parametrize("loss", [np.nan, np.inf])
def test_nonfinite_losses(data, monkeypatch, loss):
    from sklearn.linear_model import _coordinate_descent

    monkeypatch.setattr(
        _coordinate_descent, "_path_residuals", lambda *args, **kwargs: np.full(7, loss)
    )
    with pytest.raises(ValueError, match="non-finite validation losses"):
        estimator().fit(*data)


def test_metadata_routing(data):
    from sklearn import config_context

    X, y = data
    weights = np.linspace(0.5, 2, len(y))
    expected = estimator().fit(X, y, sample_weight=weights)
    with config_context(enable_metadata_routing=True):
        actual = estimator().fit(X, y, sample_weight=weights)
    assert_allclose(actual.mse_path_, expected.mse_path_)
    assert_allclose(actual.coef_, expected.coef_)


def test_overflowing_loss_moments(data, monkeypatch):
    from sklearn.linear_model import _coordinate_descent

    losses = iter([[1e200, 1]] * 3 + [[1, 1e200]] * 3)
    monkeypatch.setattr(
        _coordinate_descent,
        "_path_residuals",
        lambda *args, **kwargs: np.asarray(next(losses)),
    )
    with pytest.raises(ValueError, match="loss moments overflowed"):
        estimator(alphas=[1, 0.1], min_cv_repeats=2).fit(*data)


def test_float32(data):
    X, y = (array.astype(np.float32) for array in data)
    model = estimator(tol=1e-5).fit(X, y)
    fixed = clone(model).set_params(cv_tolerance=None).fit(X, y)
    assert_allclose(model.mse_path_, fixed.mse_path_)
    assert model.coef_.dtype == np.float32
    curves = model.mse_path_.astype(np.float64).reshape(7, 4, 3).mean(axis=2).T
    assert_allclose(
        model.cv_mcse_, (curves - curves[:, :1]).std(axis=0, ddof=1) / 2, atol=1e-14
    )
