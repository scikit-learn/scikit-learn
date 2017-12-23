"""
The :mod:`sklearn.benchmark` module includes utilities for documenting
the estimated runtime or space complexities of estimators and algorithms.
"""

# Author: Vrishank Bhardwaj <vrishank1997@gmail.com>
# License: BSD 3 clause

import time

import numpy as np

from .model_selection import train_test_split
from .utils.validation import _num_samples
from .gaussian_process import GaussianProcessRegressor
from .gaussian_process.kernels import DotProduct

__all__ = ['benchmark_estimator_cost']


def benchmark_estimator_cost(est, X, y=None, est_params=None, fit_params=None,
                             vary_n_samples=True, n_fits=5, time_budget=300,
                             profile_memory=True):
    """Profiles the cost of fitting est on samples of different size

    Parameters
    ----------
    est : estimator
    X : array-like
    y : array-like, optional
    fit_params : dict, optional
    vary_n_samples : bool, default=True
        Whether to benchmark for various random sample sizes.
    n_fits : int, default=100
        Maximum number of fits to make while benchmarking.
    time_budget : int, default=300
        Maximum number of seconds to use overall.  Current fit will
        be stopped if the budget is exceeded.
    profile_memory : bool, default=True
        Whether to include memory (or just time) profiling. Memory
        profiling will slow down fitting, and hence make fit_time
        estimates more approximate.

    Returns
    -------
    results : dict
        The following keys are each mapped to an array:

        n_samples
            The number of samples
        fit_time
            In seconds
        peak_memory
            The memory used at peak of fitting, in KiB.
        model_memory
            The memory in use at the end of fitting, minus that at the
            beginning, in KiB.

    models : dict
        keys 'peak_memory', 'model_memory' and 'fit_time' map to polynomial
        GP regressors whose input is n_samples  and whose
        outputs are each of those targets.

    errors : list of dicts
        lists the parameters that resulted in exceptions

    """
    fit_times = []
    length = _num_samples(X)
    estimator = est
    parameters = est_params
    n_samples = []
    SAMPLES = 8

    if parameters is not None:
        estimator.set_params(**parameters)

    if not vary_n_samples:
        for _ in range(0, n_fits):
            if y is None:
                start_time = time.time()
                if fit_params is not None:
                    estimator.fit(X, **fit_params)
                else:
                    estimator.fit(X)
                time_taken = time.time() - start_time
                fit_times.append(time_taken)
                n_samples.append(length)
            else:
                start_time = time.time()
                if fit_params is not None:
                    estimator.fit(X, y, **fit_params)
                else:
                    estimator.fit(X, y)
                time_taken = time.time() - start_time
                fit_times.append(time_taken)

    if vary_n_samples:
        for _ in range(0, n_fits):
            if SAMPLES < length/2:
                SAMPLES *= 2
            n_samples.append(SAMPLES)
            train_size = n_samples[_] - 1
            X_vary, X_test, y_vary, y_test = train_test_split(X, y, train_size)
            if y is None:
                if fit_params is not None:
                    start_time = time.time()
                    estimator.fit(X_vary, **fit_params)
                    time_taken = time.time() - start_time
                    fit_times.append(time_taken)
                else:
                    start_time = time.time()
                    estimator.fit(X_vary)
                    time_taken = time.time() - start_time
                    fit_times.append(time_taken)
            else:
                if fit_params is not None:
                    start_time = time.time()
                    estimator.fit(X_vary, y_vary, **fit_params)
                    time_taken = time.time() - start_time
                    fit_times.append(time_taken)
                else:
                    start_time = time.time()
                    estimator.fit(X_vary, y_vary)
                    time_taken = time.time() - start_time
                    fit_times.append(time_taken)

    models = {}

    kernel = DotProduct()
    n_samples = np.array(n_samples).reshape(-1, 1)
    '''
    -----------
    TODO
    -----------

    if profile_memory:
        fit_time = GaussianProcessRegressor()
        models["peak_memory"] = peak_memory.fit(n_samples, peak_memory)

        fit_time = GaussianProcessRegressor()
        models["model_memory"] = model_memory.fit(n_samples, model_memory)

    '''
    fit_time = GaussianProcessRegressor(kernel, alpha=0.1)
    models["fit_time"] = fit_time.fit(n_samples.reshape(-1, 1), fit_times)

    results = {}
    results["n_samples"] = n_samples
    results["fit_time"] = fit_times

    return results, models
