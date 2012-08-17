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

import numpy as np
import matplotlib.pyplot as plt

from sklearn.datasets import load_digits
from sklearn.grid_search import GridSearchCV
from sklearn.tree import DecisionTreeClassifier

iris = load_digits()
X, y = iris.data, iris.target

param_grid = {'max_depth': np.arange(1, 10, 2), 'min_samples_leaf': [1, 5, 10],
              'min_samples_split': [1, 5, 10],
              'max_features': [1, 10, 30, 40, 64]}

grid_search = GridSearchCV(DecisionTreeClassifier(), param_grid=param_grid,
                            cv=3)
grid_search.fit(X, y)

results = grid_search.scores_

fig, axes = plt.subplots(2, 2)
axes = axes.ravel()

for ax, param in zip(axes, results.params):
    ax.errorbar(results.values[param], results.accumulated_mean(param, 'max'),
            yerr=results.accumulated_std(param, 'max'))
    ax.set_title(param)
plt.show()
