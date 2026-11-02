"""Experimental adaptive repeated CV for a fixed Lasso path.

This benchmark is not a public estimator API. It supports dense, finite,
single-output, unweighted regression with an intercept. Each observation of
the stopping statistic is a complete independently randomized K-fold partition.
Fold errors are averaged equally, as in LassoCV (including unequal fold sizes).

Ritzwoller and Romano (2026), https://doi.org/10.1257/aer.20231651,
motivate stopping on the MCSE of mean relative risks. We use the largest alpha
as a fixed reference and a Bonferroni allocation across the remaining alphas.
The nominal error is an asymptotic approximation, not a finite-sample bound.
Neither the criterion nor a small risk gap guarantees a stable selected model.

Run with ``python benchmarks/bench_lasso_reproducibility.py --pairs 1000``.
The fixed training and test data are reused across independent algorithm runs.
"""

import argparse
import json
from numbers import Integral, Real
from time import perf_counter

import numpy as np
from scipy.stats import norm
from threadpoolctl import threadpool_limits

from sklearn.linear_model import Lasso, lasso_path
from sklearn.model_selection import KFold
from sklearn.utils import check_random_state
from sklearn.utils.parallel import Parallel, delayed
from sklearn.utils.validation import check_X_y


def _loss_curve(X, y, alphas, folds, *, tol, max_iter):
    """Compute one complete CV curve, without a full-data refit."""
    losses = []
    for train, test in folds:
        X_mean, y_mean = X[train].mean(axis=0), y[train].mean()
        _, coefs, _ = lasso_path(
            np.asfortranarray(X[train] - X_mean),
            y[train] - y_mean,
            alphas=alphas,
            tol=tol,
            max_iter=max_iter,
        )
        residual = y[test, None] - y_mean - (X[test] - X_mean) @ coefs
        losses.append(np.mean(residual**2, axis=0))
    curve = np.mean(losses, axis=0)
    if not np.isfinite(curve).all():
        raise ValueError("The loss curve contains non-finite values.")
    return curve


def adaptive_lasso_cv(
    X,
    y,
    *,
    alphas,
    tolerance,
    error_rate=0.05,
    n_splits=5,
    min_repeats=20,
    max_repeats=1000,
    random_state=None,
    tol=1e-8,
    max_iter=10000,
):
    """Return a fitted Lasso and a separate relative-risk precision report.

    ``tolerance`` is in squared-response units and must be chosen in advance.
    ``alphas`` must be a fixed grid of at least two distinct positive values;
    it is sorted descending, so ties favor stronger regularization. Scaling,
    if needed, must be fitted within training folds; this prototype only centers.

    The report retains per-partition loss curves, MCSEs of relative risks,
    repetition count, computational cost, and budget status. ``near_best_alphas``
    identifies observed mean risks within twice the requested tolerance of the
    minimum. It is a descriptive risk-gap diagnostic, not a confidence set.
    ``partition_winner_counts`` describes single-partition winners, not the
    probability of agreement between two independent adaptive runs.
    """
    X_original = X
    X, y = check_X_y(X, y, dtype=np.float64, y_numeric=True)
    alphas = np.asarray(alphas)
    if (
        alphas.ndim != 1
        or alphas.size < 2
        or alphas.dtype.kind not in "fiu"
        or not np.isfinite(alphas).all()
        or np.any(alphas <= 0)
        or np.unique(alphas).size != alphas.size
    ):
        raise ValueError(
            "alphas must contain at least two distinct finite positive values"
        )
    alphas = np.sort(alphas.astype(float))[::-1]
    for name, value in (
        ("tolerance", tolerance),
        ("error_rate", error_rate),
        ("tol", tol),
    ):
        if (
            isinstance(value, (bool, np.bool_))
            or not isinstance(value, Real)
            or not np.isfinite(value)
            or value <= 0
        ):
            raise ValueError(f"{name} must be a finite positive real number")
    if error_rate >= 1:
        raise ValueError("error_rate must be less than one")
    for name, value in (
        ("min_repeats", min_repeats),
        ("max_repeats", max_repeats),
        ("n_splits", n_splits),
        ("max_iter", max_iter),
    ):
        if isinstance(value, (bool, np.bool_)) or not isinstance(value, Integral):
            raise ValueError(f"{name} must be an integer")
    if not 2 <= min_repeats <= max_repeats:
        raise ValueError("Require 2 <= min_repeats <= max_repeats")
    if not 2 <= n_splits <= len(y) or max_iter < 1:
        raise ValueError("Require 2 <= n_splits <= n_samples and max_iter >= 1")
    threshold = tolerance / (
        np.sqrt(2) * norm.isf(error_rate / (2 * (alphas.size - 1)))
    )
    if not np.isfinite(threshold) or threshold <= 0:
        raise ValueError("The target is outside floating-point precision")
    splitter = KFold(
        n_splits, shuffle=True, random_state=check_random_state(random_state)
    )
    curves = []
    relative_mean = np.zeros(alphas.size)
    relative_m2 = np.zeros(alphas.size)
    criterion_met = False
    for count in range(1, max_repeats + 1):
        curve = _loss_curve(X, y, alphas, splitter.split(X), tol=tol, max_iter=max_iter)
        curves.append(curve)
        relative = curve - curve[0]
        delta = relative - relative_mean
        relative_mean += delta / count
        relative_m2 += delta * (relative - relative_mean)
        if not np.isfinite(relative_mean).all() or not np.isfinite(relative_m2).all():
            raise ValueError("Relative-loss moments overflowed; rescale the response")
        if count >= min_repeats:
            mcse = np.sqrt(np.maximum(relative_m2, 0) / (count * (count - 1)))
            criterion_met = bool(np.all(mcse[1:] <= threshold))
            if criterion_met:
                break
    curves = np.asarray(curves)
    with np.errstate(over="ignore", invalid="ignore"):
        mean_loss = curves.mean(axis=0)
    if not np.isfinite(mean_loss).all():
        raise ValueError("Mean validation loss overflowed; rescale the response")
    best = np.argmin(mean_loss)
    fitted = Lasso(alpha=alphas[best], tol=tol, max_iter=max_iter).fit(X_original, y)
    return fitted, {
        "status": "target_reached" if criterion_met else "budget_exhausted",
        "n_repeats": count,
        "n_path_fits": count * n_splits,
        "n_full_fits": 1,
        "alphas": alphas,
        "loss_curves": curves,
        "mean_loss": mean_loss,
        "relative_mean": relative_mean,
        "mcse": mcse,
        "threshold": threshold,
        "tolerance": tolerance,
        "error_rate": error_rate,
        "near_best_alphas": alphas[mean_loss - mean_loss[best] <= 2 * tolerance],
        "partition_winner_counts": np.bincount(
            np.argmin(curves, axis=1), minlength=alphas.size
        ),
    }


def _fixed_cv(X, y, alphas, repeats, seed, *, shuffled=True):
    rng = np.random.RandomState(seed)
    curves = [
        _loss_curve(
            X,
            y,
            alphas,
            KFold(5, shuffle=shuffled, random_state=rng if shuffled else None).split(X),
            tol=1e-8,
            max_iter=10000,
        )
        for _ in range(repeats)
    ]
    mean = np.mean(curves, axis=0)
    fitted = Lasso(alpha=alphas[np.argmin(mean)], tol=1e-8, max_iter=10000).fit(X, y)
    return fitted, {
        "relative_mean": mean - mean[0],
        "n_repeats": repeats,
        "n_path_fits": 5 * repeats,
    }


def _alo_cv(X, y, alphas):
    """Active-set ALO, with an unpenalized intercept and fixed sklearn alpha.

    Uses Wang et al. (2018), Section 4.2:
    https://proceedings.mlr.press/v80/wang18m.html. The full-data path uses
    alpha * (n - 1) / n so deletion holds the mean-loss penalty alpha fixed.
    Active-set changes invalidate exactness. A leverage near one is rejected.
    """
    n = len(y)
    X_mean, y_mean = X.mean(axis=0), y.mean()
    _, coefs, _ = lasso_path(
        np.asfortranarray(X - X_mean),
        y - y_mean,
        alphas=alphas * (n - 1) / n,
        tol=1e-8,
        max_iter=10000,
    )
    risks = []
    for coef in coefs.T:
        active_design = np.column_stack((np.ones(n), X[:, coef != 0]))
        leverage = np.sum(active_design * np.linalg.pinv(active_design).T, axis=1)
        if np.any(1 - leverage <= 1e-10):
            raise ValueError("ALO is undefined for leverage near one")
        residual = y - y_mean - (X - X_mean) @ coef
        risks.append(np.mean((residual / (1 - leverage)) ** 2))
    risks = np.asarray(risks)
    fitted = Lasso(alpha=alphas[np.argmin(risks)], tol=1e-8, max_iter=10000).fit(X, y)
    return fitted, {
        "relative_mean": risks - risks[0],
        "mean_loss": risks,
        "n_repeats": 0,
        "n_path_fits": 1,
    }


def _run_pair(run, pair, X_test, y_test, tolerance):
    with threadpool_limits(limits=1):
        (first, a), (second, b) = run(2 * pair), run(2 * pair + 1)
        return {
            "exceeded": np.max(np.abs(a["relative_mean"] - b["relative_mean"]))
            > tolerance,
            "alpha_agreement": first.alpha == second.alpha,
            "support_agreement": np.array_equal(first.coef_ != 0, second.coef_ != 0),
            "risk": [
                float(np.mean((y_test - f.predict(X_test)) ** 2))
                for f in (first, second)
            ],
            "repeats": [r["n_repeats"] for r in (a, b)],
            "path_fits": [r["n_path_fits"] for r in (a, b)],
            "exhausted": [r.get("status") == "budget_exhausted" for r in (a, b)],
        }


def benchmark(pairs=1000, tolerance=0.02, max_repeats=1000, n_jobs=1):
    """Measure pairwise reproducibility, prediction, selection, and work."""
    rng = np.random.RandomState(1234)
    X = rng.normal(size=(120, 20))
    X_test = rng.normal(size=(5000, 20))
    coef = np.r_[0.5, -0.4, 0.3, np.zeros(17)]
    y = X @ coef + rng.normal(size=len(X))
    y_test = X_test @ coef + rng.normal(size=len(X_test))
    alphas = np.geomspace(1, 0.01, 20)
    methods = {
        "single_randomized": lambda seed: _fixed_cv(X, y, alphas, 1, seed),
        "fixed_10": lambda seed: _fixed_cv(X, y, alphas, 10, seed),
        "fixed_50": lambda seed: _fixed_cv(X, y, alphas, 50, seed),
        "fixed_100": lambda seed: _fixed_cv(X, y, alphas, 100, seed),
        "fixed_500": lambda seed: _fixed_cv(X, y, alphas, 500, seed),
        "adaptive": lambda seed: adaptive_lasso_cv(
            X,
            y,
            alphas=alphas,
            tolerance=tolerance,
            max_repeats=max_repeats,
            random_state=seed,
        ),
        "deterministic": lambda seed: _fixed_cv(X, y, alphas, 1, seed, shuffled=False),
        "approximate_loo": lambda seed: _alo_cv(X, y, alphas),
    }
    results = {}
    for name, run in methods.items():
        start = perf_counter()
        rows = Parallel(n_jobs=n_jobs)(
            delayed(_run_pair)(run, pair, X_test, y_test, tolerance)
            for pair in range(pairs)
        )
        probability = float(np.mean([r["exceeded"] for r in rows]))
        z = norm.isf(0.025)
        center = (probability + z**2 / (2 * pairs)) / (1 + z**2 / pairs)
        half = (
            z
            * np.sqrt(probability * (1 - probability) / pairs + z**2 / (4 * pairs**2))
            / (1 + z**2 / pairs)
        )
        results[name] = {
            "pairs": pairs,
            "risk_curve_exceedance": probability,
            "exceedance_mcse": float(np.sqrt(probability * (1 - probability) / pairs)),
            "exceedance_interval_95": [float(center - half), float(center + half)],
            "alpha_agreement": float(np.mean([r["alpha_agreement"] for r in rows])),
            "support_agreement": float(np.mean([r["support_agreement"] for r in rows])),
            "test_mse": float(np.mean([r["risk"] for r in rows])),
            "mean_repeats": float(np.mean([r["repeats"] for r in rows])),
            "total_path_fits": int(np.sum([r["path_fits"] for r in rows])),
            "budget_exhaustion_rate": float(np.mean([r["exhausted"] for r in rows])),
            "elapsed_seconds": perf_counter() - start,
        }
    return {
        "seed": 1234,
        "n_jobs": n_jobs,
        "tolerance": tolerance,
        "max_repeats": max_repeats,
        "results": results,
    }


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--pairs", type=int, default=1000)
    parser.add_argument("--tolerance", type=float, default=0.02)
    parser.add_argument("--max-repeats", type=int, default=1000)
    parser.add_argument("--n-jobs", type=int, default=1)
    args = parser.parse_args()
    if args.pairs < 1:
        parser.error("--pairs must be positive")
    with threadpool_limits(limits=1):
        print(
            json.dumps(
                benchmark(args.pairs, args.tolerance, args.max_repeats, args.n_jobs),
                indent=2,
            )
        )
