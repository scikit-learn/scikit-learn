"""
==================================================================
Plotting Performance of different hyperparameter selection methods
==================================================================

In this plot you can see the validation scores of a Ridge regression model
along steps of different hyperparameter selection methods.
"""
print(__doc__)

import matplotlib.pyplot as plt

import numpy as np
from scipy import stats

from sklearn.model_selection import SequentialSearchCV
from sklearn.model_selection import GridSearchCV
from sklearn.model_selection import RandomizedSearchCV
from sklearn.linear_model import Ridge
from sklearn.datasets import make_regression


# Loading the Digits dataset
X, y = make_regression(
    n_samples=50, n_features=500, noise=1., random_state=0)

# build a classifier
clf = Ridge()

# specify parameters and distributions to sample from
expon = stats.expon(scale=20)
expon.random_state = 0  # for reproductible results
params = {'alpha': expon}

# run randomized search
gp_search = SequentialSearchCV(
  clf, params, verbose=1, n_iter=20)
gp_search.fit(X, y)

# retrieve the maximum score at each iteration
gp_cum_scores = [np.max(
      [gp_search.grid_scores_[j].mean_validation_score for j in range(i)])
      for i in range(1, len(gp_search.grid_scores_))]

rdm_search = RandomizedSearchCV(
  clf, params, random_state=0, n_iter=20)
rdm_search.fit(X, y)

rdm_cum_scores = [np.max(
      [rdm_search.grid_scores_[j].mean_validation_score for j in range(i)])
      for i in range(1, len(rdm_search.grid_scores_))]


params = {'alpha': np.logspace(-3, 3, 20)}

grid_search = GridSearchCV(clf, params)
grid_search.fit(X, y)

grid_cum_scores = [np.max(
      [grid_search.grid_scores_[j].mean_validation_score for j in range(i)])
      for i in range(1, len(grid_search.grid_scores_))]

plt.plot(gp_cum_scores, label='GP', lw=3)
plt.plot(rdm_cum_scores, label='Random Search', lw=3)
plt.plot(grid_cum_scores, label='Grid Search', lw=3)

plt.legend(loc='lower right')
plt.ylabel('Score (higher is better)')
plt.xlabel('Number of model evaluations')
plt.show()
