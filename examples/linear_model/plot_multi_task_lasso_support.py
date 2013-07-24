#!/usr/bin/env python
"""
=============================================
Joint feature selection with multi-task Lasso
=============================================

The multi-task lasso allows to fit multiple regression problems
jointly enforcing the selected features to be the same across
tasks. This example simulates sequential measurements, each task
is a time instant, and the relevant features vary in amplitude
over time while being the same. The multi-task lasso imposes that
features that are selected at one time point are select for all time
point. This makes feature selection by the Lasso more stable.

"""
print(__doc__)

# Author: Alexandre Gramfort <alexandre.gramfort@inria.fr>
# License: BSD 3 clause

import pylab as pl
import numpy as np

from sklearn.linear_model import MultiTaskLasso, Lasso

rng = np.random.RandomState(42)

# Generate some 2D coefficients with sine waves with random frequency and phase
n_samples, n_features, n_tasks = 100, 30, 40
n_relevant_features = 5
coef = np.zeros((n_tasks, n_features))
times = np.linspace(0, 2 * np.pi, n_tasks)
for k in range(n_relevant_features):
    coef[:, k] = np.sin((1. + rng.randn(1)) * times + 3 * rng.randn(1))

X = rng.randn(n_samples, n_features)
Y = np.dot(X, coef.T) + rng.randn(n_samples, n_tasks)

coef_lasso_ = np.array([Lasso(alpha=0.5).fit(X, y).coef_ for y in Y.T])
coef_multi_task_lasso_ = MultiTaskLasso(alpha=1.).fit(X, Y).coef_

###############################################################################
# Plot support and time series
fig = pl.figure(figsize=(8, 5))
pl.subplot(1, 2, 1)
pl.spy(coef_lasso_)
pl.xlabel('Feature')
pl.ylabel('Time (or Task)')
pl.text(10, 5, 'Lasso')
pl.subplot(1, 2, 2)
pl.spy(coef_multi_task_lasso_)
pl.xlabel('Feature')
pl.ylabel('Time (or Task)')
pl.text(10, 5, 'MultiTaskLasso')
fig.suptitle('Coefficient non-zero location')

feature_to_plot = 0
pl.figure()
pl.plot(coef[:, feature_to_plot], 'k', label='Ground truth')
pl.plot(coef_lasso_[:, feature_to_plot], 'g', label='Lasso')
pl.plot(coef_multi_task_lasso_[:, feature_to_plot],
        'r', label='MultiTaskLasso')
pl.legend(loc='upper center')
pl.axis('tight')
pl.ylim([-1.1, 1.1])
pl.show()
