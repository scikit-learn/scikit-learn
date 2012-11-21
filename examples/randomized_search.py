"""
========================================================================
Comparing randomized search and grid search for hyperparmeter estimation
========================================================================

Compare randomized search and grid search for optimizing hyper parameters of a
random forest.
All parameters that influence the learning are searched simultaneously
(except for the number of estimators, which poses a time / quality tradeoff).

The randomized search and the grid search explore exactly the same space of
parameters.  The result in parameter settings is quite similar, while the run
time for randomized search is drastically lower.

The performance is slightly worse for the randomized search, though this
is most likely a noise effect and would not carry over to a hold-out test set.

Note that in practice, one would not search over this many different parameters
simulatneously using grid search, but pick only the ones deemed most important.
"""
print __doc__

from time import time
from scipy.stats.distributions import randint

from sklearn.grid_search import GridSearchCV, RandomizedSearchCV
from sklearn.datasets import load_iris
from sklearn.ensemble import RandomForestClassifier

# get some data
iris = load_iris()
X, y = iris.data, iris.target

# build a classifier
clf = RandomForestClassifier(n_estimators=20)

# specify parameters and distributions to sample from
param_dist = {"max_depth": randint(1, 5), "max_features": randint(1, 4),
    "min_samples_split": randint(1, 5), "min_samples_leaf": randint(1, 5),
    "bootstrap": [True, False], "criterion": ["gini",  "entropy"]}

# run randomized search
random_search = RandomizedSearchCV(clf, param_distributions=param_dist,
        n_iter=20)

start = time()
random_search.fit(X, y)
print("RandomizedSearchCV took %.2f seconds for 20 iterations."
        % (time() - start))
print("Best score: %f" % random_search.best_score_)
print("Best parameters: %s" % repr(random_search.best_params_))

# use a full grid over all parameters
param_grid = {"max_depth": range(1, 5), "max_features": range(1, 4),
    "min_samples_split": range(1, 5), "min_samples_leaf": range(1, 5),
    "bootstrap": [True, False], "criterion": ["gini",  "entropy"]}

# run grid search
grid_search = GridSearchCV(clf, param_grid=param_grid)
grid_search.fit(X, y)

print("GridSearchCV took %.2f seconds." % (time() - start))
print("Best score: %f" % grid_search.best_score_)
print("Best parameters: %s" % repr(grid_search.best_params_))
