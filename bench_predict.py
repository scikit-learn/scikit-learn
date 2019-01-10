"""
Compare prediction time with pygbm.

run with
export NUMBA_NUM_THREADS=1 && make in && python bench_predict.py
"""

from time import time
from collections import defaultdict

import pygbm
import numpy as np
import matplotlib.pyplot as plt

from sklearn.datasets import make_regression, make_classification
from sklearn.ensemble import GradientBoostingRegressor
from sklearn.ensemble import GradientBoostingClassifier
from sklearn.ensemble import GBMRegressor
from sklearn.ensemble import GBMClassifier

classif = False
n_classes = 3
max_pow = 7
n_samples = int(10**max_pow)
max_iter = 20
n_features = 5

if classif:
    X, y = make_classification(n_samples=n_samples, n_features=n_features,
                               random_state=0, n_classes=n_classes,
                               n_clusters_per_class=1)
    GBM = GBMClassifier
    GBDT = GradientBoostingClassifier
    PYGBM_GBM = pygbm.GradientBoostingClassifier
else:
    X, y = make_regression(n_samples=n_samples, n_features=n_features,
                           random_state=0)
    GBM = GBMRegressor
    GBDT = GradientBoostingRegressor
    PYGBM_GBM = pygbm.GradientBoostingRegressor


sklearn_est = GBM(
    max_iter=max_iter,
    scoring=None,  # no early stopping
    validation_split=None,
    n_iter_no_change=None,
    random_state=0,
    verbose=False)

pygbm_est = PYGBM_GBM(
    max_iter=max_iter,
    scoring=None,  # no early stopping
    validation_split=None,
    random_state=0,
    verbose=False)
print("compiling pygbm code, and fit estimators")
pygbm_est.fit(X[:1000], y[:1000])
pygbm_est.predict(X[:1000])
sklearn_est.fit(X[:1000], y[:1000])
print("done")

n_samples_list = [10**x for x in range(2, max_pow + 1)]
n_exp = 3

predict_durations = defaultdict(lambda: defaultdict(list))

for n_samples in n_samples_list:
    for exp in range(n_exp):

        tic = time()
        sklearn_est.predict(X[:n_samples])
        predict_duration = time() - tic
        print(f'sklearn_est predict_duration: {predict_duration:.3f}s')

        predict_durations['sklearn'][n_samples].append(predict_duration)

        tic = time()
        pygbm_est.predict(X[:n_samples])
        predict_duration = time() - tic
        print(f'pygbm_est predict_duration: {predict_duration:.3f}s\n')
        predict_durations['pygbm'][n_samples].append(predict_duration)


fig, ax = plt.subplots(1)

for implem in ('sklearn', 'pygbm'):
    avgs = [np.mean(predict_durations[implem][n_samples])
            for n_samples in n_samples_list]
    stds = [np.std(predict_durations[implem][n_samples])
            for n_samples in n_samples_list]
    ax.errorbar(n_samples_list, avgs, yerr=stds, label=implem)
ax.set_xscale('log')
ax.legend(loc='best')

fig.suptitle(f'Avg prediction time over {n_exp} runs\nfor different sample sizes')
plt.show()
