"""
=====================================================
Visualizing results of high dimensional grid searches
=====================================================

Often one is faced with combining feature extraction, feature selection
and classification into a complex pipeline.
Each individual step usually has many tunable parameters.  Finding the
important parameters for a given task and picking robust settings is often
hard.

This example show how to visualize results of a grid search with
many interacting parameters.
The ``DecisionTreeClassifier`` is a good model for a complex pipeline as there
are many parameters to tweak, but only few have significant influence.
"""
print __doc__

import matplotlib.pyplot as plt

from sklearn.datasets import make_classification
from sklearn.grid_search import GridSearchCV
from sklearn.tree import DecisionTreeClassifier

X, y = make_classification(n_samples=100, n_features=10)

param_grid = {'max_depth': range(1, 8), 'min_samples_leaf': [1, 2, 3, 4, 5],
        'max_features': [1, 3, 5, 8, 10]}

grid_search = GridSearchCV(DecisionTreeClassifier(), param_grid=param_grid,
                            cv=5)
grid_search.fit(X, y)

results = grid_search.scores_

fig, axes = plt.subplots(1, 3)
axes = axes.ravel()

for ax, param in zip(axes, results.params):
    means, errors = results.accumulated(param, 'mean')
    ax.errorbar(results.values[param], means, yerr=errors)
    ax.set_title(param)
plt.show()
