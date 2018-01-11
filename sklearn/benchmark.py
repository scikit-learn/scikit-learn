"""
The :mod:`sklearn.benchmark` module includes utilities for estimating
runtimes or space complexities of estimators and algorithms.
"""

# Author: Vrishank Bhardwaj <vrishank1997@gmail.com>
# License: BSD 3 clause

import time

import numpy as np

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

    Notes
    -----
    Memory estimates are only for the current process, memory will not be
    estimated correctly if `n_jobs` with multiprocessing is used in the
    estimator.

    Peak memory measurements are approximate, as they are based on polling
    current memory usage during execution.

    """
    fit_times = []
    length = _num_samples(X)
    parameters = est_params
    n_samples = []
    model_memory = []
    peak_memory = []
    models = {}
    results = {}
    errors = {}
    errors["n_samples"] = []
    SAMPLES = 8

    if profile_memory:
        try:
            from memory_profiler import memory_usage
        except ImportError:
            raise ImportError("To run the profiler with `profile_memory=True,`"
                              " you need to install 'memory_profiler`. Please "
                              "run `pip install memory_profiler`.")

    if parameters is not None:
        est.set_params(**parameters)
    if fit_params is None:
        fit_params = {}

    for _ in range(0, n_fits):
        # if vary_n_samples=True, number of samples starts at 8 and doubles
        # for every fit.
        if vary_n_samples:
            if 8*2**_ < length:
                SAMPLES = 8*2**_
            else:
                SAMPLES = length

            n_samples.append(SAMPLES)
            train_size = n_samples[_] - 1
            X_fit = X[:train_size]
            y_fit = y[:train_size]
        # if vary_n_samples=False, the estimator is fit on the entire data for
        # every fit.
        else:
            n_samples.append(length)
            X_fit = X[:length]
            y_fit = y[:length]

        start_time = time.time()

        try:
            if profile_memory:
                args = (est, X_fit, y_fit, fit_params)
                mem = memory_usage((est_arg_handler, args), interval=.01)
                model_memory.append(mem[-1]-mem[0])
                peak_memory.append(np.max(mem))
            else:
                est_arg_handler(est, X_fit, y_fit, fit_params)
        except:
            errors["n_samples"].append(SAMPLES)

        time_taken = time.time() - start_time
        fit_times.append(time_taken)

    # Generate models

    kernel = DotProduct()
    n_samples = np.array(n_samples).reshape(-1, 1)

    if profile_memory:
        model_peak_mem = GaussianProcessRegressor()
        models["peak_memory"] = model_peak_mem.fit(n_samples, peak_memory)
        results["peak_memory"] = peak_memory

        model_mem = GaussianProcessRegressor()
        models["model_memory"] = model_mem.fit(n_samples, model_memory)
        results["model_memory"] = model_memory

    fit_time = GaussianProcessRegressor(kernel, alpha=0.1)
    models["fit_time"] = fit_time.fit(n_samples.reshape(-1, 1), fit_times)

    results["n_samples"] = n_samples
    results["fit_time"] = fit_times

    return results, models, errors


def est_arg_handler(estimator, X, y=None, fit_params=None):
    estimator.fit(X, y, **(fit_params))
    return
