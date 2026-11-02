import importlib.util
from pathlib import Path

import numpy as np
import pytest
from numpy.testing import assert_allclose, assert_array_equal
from scipy.stats import norm

from sklearn.linear_model import Lasso, LassoCV
from sklearn.model_selection import KFold, LeaveOneOut, RepeatedKFold
from sklearn.utils.validation import check_is_fitted

spec = importlib.util.spec_from_file_location(
    "lasso_reproducibility",
    Path(__file__).parents[1] / "bench_lasso_reproducibility.py",
)
benchmark = importlib.util.module_from_spec(spec)
spec.loader.exec_module(benchmark)


@pytest.fixture
def data():
    rng = np.random.RandomState(12)
    X = rng.normal(size=(53, 6)) + 3
    y = 2 + X[:, 0] + rng.normal(size=len(X))
    return X, y


def test_loss_curve_matches_lassocv(data):
    X, y = data
    alphas = np.geomspace(1, 0.01, 7)
    folds = list(KFold(5, shuffle=True, random_state=20).split(X))
    original_X, original_y = X.copy(), y.copy()
    actual = benchmark._loss_curve(X, y, alphas, folds, tol=1e-10, max_iter=10000)
    expected = LassoCV(alphas=alphas, cv=folds, tol=1e-10, max_iter=10000).fit(X, y)
    assert_allclose(actual, expected.mse_path_.mean(axis=1), rtol=1e-8)
    assert_array_equal(X, original_X)
    assert_array_equal(y, original_y)


def test_loss_curve_matches_individual_fits(data):
    X, y = data
    alphas = np.geomspace(1, 0.01, 7)
    folds = list(KFold(5, shuffle=True, random_state=20).split(X))
    expected = []
    for train, test in folds:
        losses = []
        for alpha in alphas:
            fitted = Lasso(alpha=alpha, tol=1e-12, max_iter=10000).fit(
                X[train], y[train]
            )
            losses.append(np.mean((y[test] - fitted.predict(X[test])) ** 2))
        expected.append(losses)
    actual = benchmark._loss_curve(X, y, alphas, folds, tol=1e-12, max_iter=10000)
    assert_allclose(actual, np.mean(expected, axis=0), rtol=1e-10)


@pytest.mark.parametrize(
    "exceeded, interval",
    [(False, [0, 0.2775327998628892]), (True, [0.7224672001371107, 1])],
)
def test_benchmark_exceedance_interval(exceeded, interval, monkeypatch):
    def run_pair(*args):
        return {
            "exceeded": exceeded,
            "alpha_agreement": True,
            "support_agreement": True,
            "risk": [1, 1],
            "repeats": [20, 20],
            "path_fits": [100, 100],
            "exhausted": [False, False],
        }

    monkeypatch.setattr(benchmark, "_run_pair", run_pair)
    report = benchmark.benchmark(pairs=10)
    for result in report["results"].values():
        assert result["risk_curve_exceedance"] == float(exceeded)
        assert_allclose(result["exceedance_interval_95"], interval, atol=1e-15)


def test_fixed_repeated_cv_equivalence(data):
    X, y = data
    alphas = np.geomspace(1, 0.01, 7)
    fitted, report = benchmark.adaptive_lasso_cv(
        X,
        y,
        alphas=alphas,
        tolerance=1e-14,
        min_repeats=4,
        max_repeats=4,
        random_state=10,
    )
    expected = LassoCV(
        alphas=alphas,
        cv=RepeatedKFold(n_splits=5, n_repeats=4, random_state=10),
        tol=1e-8,
        max_iter=10000,
    ).fit(X, y)
    assert_allclose(report["mean_loss"], expected.mse_path_.mean(axis=1), rtol=1e-7)
    assert fitted.alpha == expected.alpha_
    assert_allclose(fitted.coef_, expected.coef_, rtol=1e-6, atol=1e-8)
    assert report["status"] == "budget_exhausted"
    assert report["n_repeats"] == 4


def test_mcse_and_bonferroni_oracle(data):
    _, report = benchmark.adaptive_lasso_cv(
        *data,
        alphas=[0.01, 0.1, 1],
        tolerance=0.01,
        min_repeats=3,
        max_repeats=6,
        random_state=12,
    )
    relative = report["loss_curves"] - report["loss_curves"][:, :1]
    assert_allclose(report["relative_mean"], relative.mean(axis=0))
    assert_allclose(
        report["mcse"], relative.std(axis=0, ddof=1) / np.sqrt(len(relative))
    )
    assert_allclose(report["threshold"], 0.01 / (np.sqrt(2) * norm.isf(0.05 / 4)))
    assert report["mcse"][0] == 0


def test_work_count_and_reproducibility(data, monkeypatch):
    calls, refits = [], []
    path = benchmark.lasso_path
    fit = benchmark.Lasso.fit

    def counted_path(*args, **kwargs):
        calls.append(kwargs["alphas"].copy())
        return path(*args, **kwargs)

    def counted_fit(self, *args, **kwargs):
        refits.append(self.alpha)
        return fit(self, *args, **kwargs)

    monkeypatch.setattr(benchmark, "lasso_path", counted_path)
    monkeypatch.setattr(benchmark.Lasso, "fit", counted_fit)
    _, report = benchmark.adaptive_lasso_cv(
        *data,
        alphas=[0.01, 1, 0.1],
        tolerance=100,
        min_repeats=3,
        max_repeats=10,
        random_state=12,
    )
    assert report["status"] == "target_reached"
    assert len(calls) == report["n_path_fits"] == 15
    assert len(refits) == report["n_full_fits"] == 1
    for grid in calls:
        assert_array_equal(grid, [1, 0.1, 0.01])
    _, again = benchmark.adaptive_lasso_cv(
        *data,
        alphas=[0.01, 1, 0.1],
        tolerance=100,
        min_repeats=3,
        max_repeats=10,
        random_state=12,
    )
    assert_array_equal(report["loss_curves"], again["loss_curves"])


def test_common_noise_cancels_and_mean_precedes_argmin(data, monkeypatch):
    curves = iter([np.array([10, 9, 11]), np.array([20, 19, 18])])
    monkeypatch.setattr(benchmark, "_loss_curve", lambda *a, **kw: next(curves))
    fitted, report = benchmark.adaptive_lasso_cv(
        *data,
        alphas=[1, 0.1, 0.01],
        tolerance=100,
        min_repeats=2,
        max_repeats=2,
    )
    assert report["mcse"][1] == 0
    assert fitted.alpha == 0.1
    assert report["partition_winner_counts"].tolist() == [0, 1, 1]
    assert report["near_best_alphas"].size == 3


@pytest.mark.parametrize(
    "alphas",
    [[], [1], [1, 1], [0, 1], [-1, 1], [np.nan, 1], [np.inf, 1], [[1, 2]], ["1", "2"]],
)
def test_invalid_grid(data, alphas):
    with pytest.raises(ValueError):
        benchmark.adaptive_lasso_cv(*data, alphas=alphas, tolerance=0.01)


@pytest.mark.parametrize(
    "kwargs",
    [
        {"tolerance": 0},
        {"tolerance": np.nan},
        {"tolerance": True},
        {"error_rate": 1},
        {"error_rate": 0},
        {"error_rate": np.inf},
        {"min_repeats": 1},
        {"min_repeats": 2.1},
        {"max_repeats": 2},
        {"n_splits": 100},
        {"n_splits": True},
        {"max_iter": 0},
        {"tol": -1},
    ],
)
def test_invalid_settings(data, kwargs):
    settings = {"alphas": [1, 0.1], "tolerance": 0.01, **kwargs}
    with pytest.raises(ValueError):
        benchmark.adaptive_lasso_cv(*data, **settings)


def test_zero_variance(data, monkeypatch):
    monkeypatch.setattr(benchmark, "_loss_curve", lambda *a, **kw: np.array([1.0, 1.0]))
    fitted, report = benchmark.adaptive_lasso_cv(
        *data,
        alphas=[0.1, 1],
        tolerance=0.001,
        min_repeats=3,
        max_repeats=3,
    )
    assert fitted.alpha == 1
    assert report["status"] == "target_reached"
    assert_array_equal(report["mcse"], [0, 0])


def test_alo_matches_loo_with_constant_active_set():
    rng = np.random.RandomState(18)
    X = rng.normal(size=(60, 1))
    y = 2 + 4 * X[:, 0] + rng.normal(scale=0.05, size=len(X))
    alphas = np.array([0.2, 0.1])
    _, report = benchmark._alo_cv(X, y, alphas)
    exact = []
    for alpha in alphas:
        residuals = []
        for train, test in LeaveOneOut().split(X):
            fitted = Lasso(alpha=alpha, tol=1e-12, max_iter=10000).fit(
                X[train], y[train]
            )
            assert fitted.coef_[0] > 0
            residuals.append((y[test] - fitted.predict(X[test])) ** 2)
        exact.append(np.mean(residuals))
    assert_allclose(report["mean_loss"], exact, rtol=1e-8)


def test_alo_null_model(data):
    X, y = data
    _, report = benchmark._alo_cv(X, y, np.array([1000.0, 500.0]))
    expected = np.mean(((y - y.mean()) / (1 - 1 / len(y))) ** 2)
    assert_allclose(report["mean_loss"], expected)


def test_overflow_cannot_certify_precision(data, monkeypatch):
    monkeypatch.setattr(
        benchmark, "_loss_curve", lambda *a, **kw: np.array([1e307, 1e307])
    )
    with pytest.raises(ValueError, match="Mean validation loss overflowed"):
        benchmark.adaptive_lasso_cv(
            *data,
            alphas=[1, 0.1],
            tolerance=0.01,
            min_repeats=20,
            max_repeats=20,
        )


def test_result_contract(data):
    fitted, report = benchmark.adaptive_lasso_cv(
        *data,
        alphas=[1, 0.1, 0.01],
        tolerance=1,
        min_repeats=3,
        max_repeats=3,
        random_state=12,
    )
    assert isinstance(fitted, Lasso)
    check_is_fitted(fitted)
    assert fitted.predict(data[0]).shape == data[1].shape
    assert fitted.coef_.shape == (data[0].shape[1],)
    assert isinstance(report, dict)
    assert set(report) == {
        "status",
        "n_repeats",
        "n_path_fits",
        "n_full_fits",
        "alphas",
        "loss_curves",
        "mean_loss",
        "relative_mean",
        "mcse",
        "threshold",
        "tolerance",
        "error_rate",
        "near_best_alphas",
        "partition_winner_counts",
    }
    assert report["loss_curves"].shape == (report["n_repeats"], 3)
    for name in (
        "alphas",
        "mean_loss",
        "relative_mean",
        "mcse",
        "partition_winner_counts",
    ):
        assert isinstance(report[name], np.ndarray)
        assert report[name].shape == (3,)
    assert report["partition_winner_counts"].sum() == report["n_repeats"]
    assert fitted.alpha == report["alphas"][np.argmin(report["mean_loss"])]


def test_dataframe_feature_names_are_preserved(data):
    pd = pytest.importorskip("pandas")
    X, y = data
    X = pd.DataFrame(X, columns=[f"feature_{j}" for j in range(X.shape[1])])
    fitted, _ = benchmark.adaptive_lasso_cv(
        X,
        y,
        alphas=[1, 0.1],
        tolerance=1,
        min_repeats=3,
        max_repeats=3,
        random_state=12,
    )
    assert_array_equal(fitted.feature_names_in_, X.columns)
    assert fitted.predict(X).shape == y.shape
    with pytest.raises(ValueError, match="feature names"):
        fitted.predict(X[X.columns[::-1]])
