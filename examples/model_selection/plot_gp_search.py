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
from sklearn.linear_model import Ridge, LogisticRegression
from sklearn.datasets import load_iris, make_classification
from sklearn.gaussian_process.kernels import Matern


# Make synthetic dataset where not all the features are explanatory
iris = load_iris()
rng = np.random.RandomState(0)
X, y = iris.data, iris.target

clf = LogisticRegression()

params = {'C': {'bounds': [10**-3, 10**3], 'scale': 'log'}}

# Run SequentialSearch with first 5 iterations random.
gp_search = SequentialSearchCV(
  clf, params, n_iter=20, random_state=0, n_init=3)
gp_search.fit(X, y)

# Retrieve the maximum score at each iteration
gp_cum_scores = [np.max(
      [gp_search.grid_scores_[j].mean_validation_score for j in range(i)])
      for i in range(1, len(gp_search.grid_scores_))]

# Do the same experiment with randomized search cv
params = {'C': stats.expon(scale=10**-1)}
rdm_search = RandomizedSearchCV(
  clf, params, random_state=0, n_iter=20)
rdm_search.fit(X, y)

rdm_cum_scores = [np.max(
      [rdm_search.grid_scores_[j].mean_validation_score for j in range(i)])
      for i in range(1, len(rdm_search.grid_scores_))]

# Do a standard grid search across pre-defined parameters.
params = {'C': np.logspace(-3, 3, 20)}

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
