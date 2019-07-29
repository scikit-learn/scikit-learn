"""
Comparison between grid search and successive halving
=====================================================
"""
from time import time

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

from sklearn.svm import SVC
from sklearn import datasets
from sklearn.model_selection import GridSearchCV
from sklearn.model_selection import GridHalvingSearchCV


rng = np.random.RandomState(0)
X, y = datasets.make_classification(n_samples=1000, random_state=rng)

gammas = [1e-1, 1e-2, 1e-3, 1e-4, 1e-5, 1e-6, 1e-7]
Cs = [1, 10, 100, 1e3, 1e4, 1e5]
param_grid = {'gamma': gammas, 'C': Cs}

clf = SVC(random_state=rng)
tic = time()
gsh = GridHalvingSearchCV(
    estimator=clf,
    param_grid=param_grid,
    budget_on='n_samples',  # budget is the number of samples
    max_budget='auto',  # max_budget=n_samples
    force_exhaust_budget=True,
    cv=5,
    ratio=2,
    random_state=rng)
gsh.fit(X, y)
gsh_time = time() - tic

tic = time()
gs = GridSearchCV(
    estimator=clf,
    param_grid=param_grid,
    cv=5)
gs.fit(X, y)
gs_time = time() - tic


def make_heatmap(ax, gs, show_iter=False, make_cbar=False):
    results = pd.DataFrame.from_dict(gs.cv_results_)
    results['params_str'] = results.params.apply(str)
    scores = results.groupby(['param_gamma', 'param_C']).mean_test_score.max()
    scores_matrix = scores.values.reshape(len(gammas), len(Cs))

    im = ax.imshow(scores_matrix)

    ax.set_xticks(np.arange(len(Cs)))
    ax.set_xticklabels(['{:.0E}'.format(x) for x in Cs])
    ax.set_xlabel('C', fontsize=15)

    ax.set_yticks(np.arange(len(gammas)))
    ax.set_yticklabels(['{:.0E}'.format(x) for x in gammas])
    ax.set_ylabel('gamma', fontsize=15)

    # Rotate the tick labels and set their alignment.
    plt.setp(ax.get_xticklabels(), rotation=45, ha="right",
             rotation_mode="anchor")

    if show_iter:
        iterations = results.groupby(['param_gamma', 'param_C']).iter.max()
        iterations_matrix = iterations.values.reshape(len(gammas), len(Cs))
        for i in range(len(gammas)):
            for j in range(len(Cs)):
                ax.text(j, i, iterations_matrix[i, j],
                        ha="center", va="center", color="w", fontsize=20)

    if make_cbar:
        fig.subplots_adjust(right=0.8)
        cbar_ax = fig.add_axes([0.85, 0.15, 0.05, 0.7])
        fig.colorbar(im, cax=cbar_ax)
        cbar_ax.set_ylabel('max mean_test_score', rotation=-90, va="bottom",
                           fontsize=15)


fig, axes = plt.subplots(ncols=2)
ax1, ax2 = axes

make_heatmap(ax1, gsh, show_iter=True)
make_heatmap(ax2, gs, make_cbar=True)

ax1.set_title('Successive Halving (time = {:.3f}s)'.format(gsh_time),
              fontsize=15)
ax2.set_title('GridSearch (time = {:.3f}s)'.format(gs_time), fontsize=15)

plt.show()
