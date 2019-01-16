from collections import defaultdict
from time import time

import numpy as np
import matplotlib.pyplot as plt
from sklearn.gbm.types import HISTOGRAM_DTYPE
from sklearn.gbm.types import X_DTYPE
from sklearn.gbm.types import X_BINNED_DTYPE
from sklearn.gbm.types import Y_DTYPE
from sklearn.gbm.splitting import SplittingContext
from sklearn.gbm.splitting import find_node_split
from sklearn.gbm.splitting import split_indices
from pygbm.splitting import SplittingContext as SplittingContext_pygbm
from pygbm.splitting import find_node_split as find_node_split_pygbm
from pygbm.splitting import split_indices as split_indices_pygbm

rng = np.random.RandomState(42)

n_bins = 255
n_features = 20  # Number of features has huge impact, it's weird
l2_regularization = 0.
min_hessian_to_split = 1e-3
min_samples_leaf = 1
min_gain_to_split = 0.

max_pow = 7
n_samples_list = [10**x for x in range(2, max_pow + 1)]
n_exp = 10

n_samples = 10**max_pow

X_binned_ = rng.randint(0, n_bins, size=(n_samples, n_features), dtype=np.uint8)
sample_indices_ = np.arange(n_samples, dtype=np.uint32)
all_gradients_ = rng.randn(n_samples).astype(Y_DTYPE)
all_hessians_ = rng.lognormal(size=n_samples).astype(Y_DTYPE)

def one_run(n_samples):

    X_binned = X_binned_[:n_samples]
    X_binned = np.asfortranarray(X_binned)
    sample_indices = sample_indices_[:n_samples]
    all_gradients = all_gradients_[:n_samples]
    all_hessians = all_hessians_[:n_samples]

    n_bins_per_feature = np.array([n_bins] * X_binned.shape[1], dtype=np.uint32)

    sklearn_context = SplittingContext(X_binned, n_bins,
                            n_bins_per_feature,
                            all_gradients, all_hessians,
                            l2_regularization, min_hessian_to_split,
                            min_samples_leaf, min_gain_to_split)
    all_gradients = all_gradients.astype(np.float32)
    all_hessians = all_hessians.astype(np.float32)
    pygbm_context = SplittingContext_pygbm(X_binned, n_bins,
                                           n_bins_per_feature,
                                           all_gradients, all_hessians,
                                           l2_regularization, min_hessian_to_split,
                                           min_samples_leaf, min_gain_to_split)

    sample_indices = np.arange(n_samples, dtype=np.uint32)

    histograms = np.zeros(shape=(n_features, n_bins), dtype=HISTOGRAM_DTYPE)
    split_info = find_node_split(sklearn_context, sample_indices, histograms)
    tic = time()
    _, _, _ = split_indices(sklearn_context, split_info, sample_indices)
    sklearn_duration = time() - tic

    split_info, _ = find_node_split_pygbm(pygbm_context, sample_indices)
    tic = time()
    _, _ = split_indices_pygbm(pygbm_context, split_info, sample_indices)
    pygbm_duration = time() - tic

    return sklearn_duration, pygbm_duration

one_run(100)  # compile pygbm

durations = defaultdict(lambda: defaultdict(list))

for n_samples in n_samples_list:
    for exp in range(n_exp):

        sklearn_duration, pygbm_duration = one_run(n_samples)
        print(f"sklearn fit duration = {sklearn_duration:.3f}")
        print(f"pygbm fit duration = {pygbm_duration:.3f}")
        durations['sklearn'][n_samples].append(sklearn_duration)
        durations['pygbm'][n_samples].append(pygbm_duration)

fig, ax = plt.subplots(1)

for implem in ('sklearn', 'pygbm'):
    avgs = [np.mean(durations[implem][n_samples])
            for n_samples in n_samples_list]
    stds = [np.std(durations[implem][n_samples])
            for n_samples in n_samples_list]
    ax.errorbar(n_samples_list, avgs, yerr=stds, label=implem)


ax.set_xscale('log')
ax.legend(loc='best')

fig.suptitle(f'Avg time for split_indices over {n_exp} runs\nfor different sample sizes')
plt.show()
