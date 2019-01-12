"""
Compare binning fitting and transform time with pygbm.
"""
from time import time
from collections import defaultdict

import numpy as np
import pygbm
import matplotlib.pyplot as plt
from sklearn.datasets import make_regression

from sklearn.ensemble.gbm.binning import BinMapper


n_features = 5

max_pow = 7
n_samples = int(10**max_pow)
X, y = make_regression(n_samples=n_samples, n_features=n_features,
                       random_state=0)

print("compiling pygbm")
pygbm_bm = pygbm.binning.BinMapper()
pygbm_bm.fit_transform(X[:1000])
print('done')

bm = BinMapper()

n_samples_list = [10**x for x in range(2, max_pow + 1)]
n_exp = 10

transform_durations = defaultdict(lambda: defaultdict(list))
fit_durations = defaultdict(lambda: defaultdict(list))

for n_samples in n_samples_list:
    for exp in range(n_exp):

        tic = time()
        tic = time()
        bm.fit(X[:n_samples])
        fit_duration = time() - tic
        print(f"sklearn fit duration = {fit_duration:.3f}")
        tic = time()
        bm.transform(X[:n_samples])
        transform_duration = time() - tic
        print(f"sklearn transform duration = {transform_duration:.3f}")

        fit_durations['sklearn'][n_samples].append(fit_duration)
        transform_durations['sklearn'][n_samples].append(transform_duration)

        tic = time()
        pygbm_bm.fit(X[:n_samples])
        fit_duration = time() - tic
        print(f"pygbm fit duration = {fit_duration:.3f}")
        tic = time()
        pygbm_bm.transform(X[:n_samples])
        transform_duration = time() - tic
        print(f"pygbm transform duration = {transform_duration:.3f}")
        fit_durations['pygbm'][n_samples].append(fit_duration)
        transform_durations['pygbm'][n_samples].append(transform_duration)

fig, axs = plt.subplots(2)

for implem in ('sklearn', 'pygbm'):
    avgs = [np.mean(fit_durations[implem][n_samples])
            for n_samples in n_samples_list]
    stds = [np.std(fit_durations[implem][n_samples])
            for n_samples in n_samples_list]
    axs[0].errorbar(n_samples_list, avgs, yerr=stds, label=implem)
    axs[0].set_title('Fit')

for implem in ('sklearn', 'pygbm'):
    avgs = [np.mean(transform_durations[implem][n_samples])
            for n_samples in n_samples_list]
    stds = [np.std(transform_durations[implem][n_samples])
            for n_samples in n_samples_list]
    axs[1].errorbar(n_samples_list, avgs, yerr=stds, label=implem)
    axs[1].set_title('transform')

for ax in axs:
    ax.set_xscale('log')
    ax.legend(loc='best')

fig.suptitle(f'Avg fit and transform time for binning over {n_exp} runs\nfor different sample sizes')
plt.show()
