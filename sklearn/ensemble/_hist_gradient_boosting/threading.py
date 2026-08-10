"""Threading helpers for histogram-based gradient boosting."""

# Authors: The scikit-learn developers
# SPDX-License-Identifier: BSD-3-Clause

import os
from time import perf_counter

from threadpoolctl import threadpool_limits

from sklearn.base import clone


def _estimate_time(func, *args, min_time):
    """Estimate run time in seconds.

    Calls ``func`` until ``min_time`` seconds elapsed, return mean seconds per call.
    """
    total_time = 0.0
    n_calls = 0
    while total_time < min_time:
        start = perf_counter()
        func(*args)
        total_time += perf_counter() - start
        n_calls += 1
    return total_time / n_calls


def estimate_speedup(estimator, method, X, y, min_time=2.0):
    """Estimate OpenMP speedup relative to single-threaded runtime.

    Mean runtimes are measured by repeated calls until at `min_time` seconds
    have elapsed.

    Parameters
    ----------
    estimator : estimator object
       The estimator to time.

    method : {'fit', 'predict'}
        The method to time.

    X : array-like of shape (n_samples, n_features)
        The input samples.

    y : array-like of shape (n_samples,)
        Targets.

    min_time : float, default=2.0
        Minimal time in seconds.

    Returns
    -------
    speedup : dict[int, float]
        Return speedup versus `n_threads == 1` for each `n_threads`
        up to ``os.cpu_count()``.
    """
    if method not in {"fit", "predict"}:
        raise ValueError(f"method must be 'fit' or 'predict', got {method}")

    est = clone(estimator)
    if method == "predict":
        est.fit(X, y)

    times = {}
    max_n_threads = os.cpu_count()
    for n_threads in range(1, max_n_threads + 1):
        with threadpool_limits(limits=n_threads, user_api="openmp"):
            if method == "fit":
                times[n_threads] = _estimate_time(est.fit, X, y, min_time=min_time)
            else:
                times[n_threads] = _estimate_time(est.predict, X, min_time=min_time)

    t1 = times[1]
    return {n: t1 / t for n, t in times.items()}
