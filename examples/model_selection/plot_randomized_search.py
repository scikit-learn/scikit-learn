"""
=========================================================================
Comparing randomized search and grid search for hyperparameter estimation
=========================================================================

Compare randomized search and grid search for optimizing hyperparameters of a
random forest.
All parameters that influence the learning are searched simultaneously
(except for the number of estimators, which poses a time / quality tradeoff).

The randomized search and the grid search explore exactly the same space of
parameters. The result in parameter settings is quite similar, while the run
time for randomized search is drastically lower.

The performance is slightly worse for the randomized search, though this
is most likely a noise effect and would not carry over to a held-out test set.

Note that in practice, one would not search over this many different parameters
simultaneously using grid search, but pick only the ones deemed most important.
"""
print(__doc__)

import numpy as np

from time import time
import scipy.stats as stats
from sklearn.utils.random import loguniform

from sklearn.model_selection import GridSearchCV, RandomizedSearchCV
from sklearn.datasets import load_digits
from sklearn.linear_model import SGDClassifier

# get some data
digits = load_digits()
X, y = digits.data, digits.target

# build a classifier
clf = SGDClassifier(loss='hinge', penalty='elasticnet',
                    fit_intercept=True)


# Utility function to report best scores
def report(results, n_top=3):
    for i in range(1, n_top + 1):
        candidates = np.flatnonzero(results['rank_test_score'] == i)
        for candidate in candidates:
            print("Model with rank: {0}".format(i))
            print("Mean validation score: {0:.3f} (std: {1:.3f})".format(
                  results['mean_test_score'][candidate],
                  results['std_test_score'][candidate]))
            print("Parameters: {0}".format(results['params'][candidate]))
            print("")


# specify parameters and distributions to sample from
param_dist = {'average': [True, False],
              'l1_ratio': stats.uniform(0, 1),
              'alpha': loguniform(-4, 0, base=10)}

# run randomized search
n_iter_search = 20
random_search = RandomizedSearchCV(clf, param_distributions=param_dist,
                                   n_iter=n_iter_search, cv=5, iid=False)

start = time()
random_search.fit(X, y)
print("RandomizedSearchCV took %.2f seconds for %d candidates"
      " parameter settings." % ((time() - start), n_iter_search))
report(random_search.cv_results_)

# use a full grid over all parameters
param_grid = {'average': [True, False],
              'l1_ratio': np.linspace(0, 1, num=10),
              'alpha': [10**-i for i in [-4, -3, -2, -1, 0]]}

# run grid search
grid_search = GridSearchCV(clf, param_grid=param_grid, cv=5, iid=False)
start = time()
grid_search.fit(X, y)

print("GridSearchCV took %.2f seconds for %d candidate parameter settings."
      % (time() - start, len(grid_search.cv_results_['params'])))
report(grid_search.cv_results_)
