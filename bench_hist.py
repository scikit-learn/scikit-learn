"""
Compare histogram building function with pygbm.

might be a bit unfair to cython code since we're calling the python versions
of the cpdef functions, which causes unnecessary conversions.
"""
from time import time
from collections import defaultdict

import matplotlib.pyplot as plt
import numpy as np
from joblib import Memory
from pygbm.histogram import _build_histogram_naive as pygbm_build_histogram_naive
from pygbm.histogram import _build_histogram as pygbm_build_histogram
from pygbm.histogram import _build_histogram_no_hessian as pygbm_build_histogram_no_hessian
from pygbm.histogram import _build_histogram_root as pygbm_build_histogram_root
from pygbm.histogram import _build_histogram_root_no_hessian as pygbm_build_histogram_root_no_hessian
from pygbm.histogram import _subtract_histograms as pygbm_subtract_histograms

from sklearn.gbm.histogram import _build_histogram_naive
from sklearn.gbm.histogram import _build_histogram
from sklearn.gbm.histogram import _build_histogram_no_hessian
from sklearn.gbm.histogram import _build_histogram_root
from sklearn.gbm.histogram import _build_histogram_root_no_hessian
from sklearn.gbm.histogram import _subtract_histograms
from sklearn.gbm.types import HISTOGRAM_DTYPE
from sklearn.gbm.types import X_DTYPE
from sklearn.gbm.types import X_BINNED_DTYPE
from sklearn.gbm.types import Y_DTYPE


m = Memory(location='/tmp')

@m.cache
def make_data(n_bins=256, n_samples=int(1e8), seed=42):
    rng = np.random.RandomState(seed)

    sample_indices = np.arange(n_samples, dtype=np.uint32)
    ordered_gradients = rng.randn(n_samples).astype(Y_DTYPE)
    ordered_hessians = rng.exponential(size=n_samples).astype(Y_DTYPE)
    binned_feature = rng.randint(0, n_bins, size=n_samples, dtype=X_BINNED_DTYPE)
    return sample_indices, binned_feature, ordered_gradients, ordered_hessians


n_bins = 256
print(f"Compiling pygbm...")
sample_indices, binned_feature, gradients, hessians = make_data(
    n_bins, n_samples=10)
tic = time()
a = pygbm_build_histogram_naive(n_bins, sample_indices, binned_feature, gradients, hessians)
b = pygbm_build_histogram(n_bins, sample_indices, binned_feature, gradients, hessians)
pygbm_subtract_histograms(n_bins, a, b)
pygbm_build_histogram_no_hessian(n_bins, sample_indices, binned_feature, gradients)
pygbm_build_histogram_root(n_bins, binned_feature, gradients, hessians)
pygbm_build_histogram_root_no_hessian(n_bins, binned_feature, gradients)
toc = time()
duration = toc - tic
print(f"done in {duration:.3f}s")

def one_run(sklearn_fun, pygbm_fun):
    print('-' * 10)
    print(sklearn_fun.__name__)

    if 'subtract' in sklearn_fun.__name__:
        # specal case for subtract... crappy
        a = pygbm_build_histogram(n_bins, sample_indices, binned_feature, gradients, hessians)
        b = pygbm_build_histogram(n_bins, sample_indices, binned_feature, gradients, hessians)

        args = [n_bins, a, b]
        tic = time()
        pygbm_fun(*args)
        pygbm_duration = time() - tic
        print(f"pygbm: Built in {pygbm_duration:.3f}s")

        a = a.astype(HISTOGRAM_DTYPE)
        b = b.astype(HISTOGRAM_DTYPE)
        args = [n_bins, a, b]
        tic = time()
        histogram = np.zeros(n_bins, dtype=HISTOGRAM_DTYPE)
        args.append(histogram)
        sklearn_fun(*args)
        sklearn_duration = time() - tic
        print(f"sklearn: Built in {sklearn_duration:.3f}s")

    else:
        args = [n_bins]
        if not 'root' in sklearn_fun.__name__:
            args.append(sample_indices)
        args += [binned_feature, gradients, hessians]
        if 'no_hessian' in sklearn_fun.__name__:
            args.pop()

        tic = time()
        pygbm_fun(*args)
        pygbm_duration = time() - tic
        print(f"pygbm: Built in {pygbm_duration:.3f}s")

        tic = time()
        histogram = np.zeros(n_bins, dtype=HISTOGRAM_DTYPE)
        args.append(histogram)
        sklearn_fun(*args)
        sklearn_duration = time() - tic
        print(f"sklearn: Built in {sklearn_duration:.3f}s")

    return sklearn_duration, pygbm_duration

n_exp = 10
n_samples_list = [10**x for x in range(2, 9)]


n_rows = 3
n_cols = 2
fig, axs = plt.subplots(n_rows, n_cols, sharex=True)

for i, (sklearn_fun, pygbm_fun) in enumerate((
        (_build_histogram_naive, pygbm_build_histogram_naive),
        (_build_histogram, pygbm_build_histogram),
        (_build_histogram_no_hessian, pygbm_build_histogram_no_hessian),
        (_build_histogram_root, pygbm_build_histogram_root),
        (_build_histogram_root_no_hessian, pygbm_build_histogram_root_no_hessian),
        (_subtract_histograms, pygbm_subtract_histograms))):

    row = i // n_cols
    col = i % n_cols
    ax = axs[row][col]

    durations = defaultdict(lambda: defaultdict(list))
    for n_samples in n_samples_list:
        sample_indices, binned_feature, gradients, hessians = make_data(
            n_bins, n_samples)
        for _ in range(n_exp):
            sklearn_duration, pygbm_duration = one_run(sklearn_fun, pygbm_fun)
            durations[n_samples]['sklearn'].append(sklearn_duration)
            durations[n_samples]['pygbm'].append(pygbm_duration)

    sklearn_avgs = [np.mean(durations[n_samples]['sklearn']) for n_samples in n_samples_list]
    sklearn_stds = [np.std(durations[n_samples]['sklearn']) for n_samples in n_samples_list]
    ax.errorbar(n_samples_list, sklearn_avgs, yerr=sklearn_stds, label='PR')

    pygbm_avgs = [np.mean(durations[n_samples]['pygbm']) for n_samples in n_samples_list]
    pygbm_stds = [np.std(durations[n_samples]['pygbm']) for n_samples in n_samples_list]
    ax.errorbar(n_samples_list, pygbm_avgs, yerr=pygbm_stds, label='pygbm')
    ax.set_xscale('log')
    ax.set_title(sklearn_fun.__name__)
    ax.legend()
fig.suptitle(f'Avg histogram computation time over {n_exp} runs\nfor different sample sizes')
plt.show()
