"""
====================================================
Comparison of cross validated score with overfitting
====================================================

These two plots compare the cross validated score of a the regression of
a simple function. We see that before the maximum value of 7 the regression is
far for the real function. On the other hand, for higher number of leaves we
clearly overfit.

"""
print __doc__

import numpy as np
from sklearn import tree


def plot_pruned_path(scores, with_std=True):
    """Plots the cross validated scores versus the number of leaves of trees"""
    import matplotlib.pyplot as plt
    means = np.array([np.mean(s) for s in scores])
    stds = np.array([np.std(s) for s in scores]) / np.sqrt(len(scores[1]))

    x = range(len(scores) + 1, 1, -1)

    plt.plot(x, means)
    if with_std:
        plt.plot(x, means + 2 * stds, lw=1, c='0.7')
        plt.plot(x, means - 2 * stds, lw=1, c='0.7')

    plt.xlabel('Number of leaves')
    plt.ylabel('Cross validated score')


# Create a random dataset
rng = np.random.RandomState(1)
X = np.sort(5 * rng.rand(80, 1), axis=0)
y = np.sin(X).ravel()
y[1::5] += 3 * (0.5 - rng.rand(16))


clf = tree.DecisionTreeRegressor(max_depth=20)
scores = tree.prune_path(clf, X, y, max_n_leaves=20,
                                    n_iterations=100, random_state=0)
plot_pruned_path(scores)

clf = clf.fit(X, y)
X_test = np.arange(0.0, 5.0, 0.01)[:, np.newaxis]

#Prepare the different pruned level
clf = clf.prune(15)
y_15 = clf.predict(X_test)

clf = clf.prune(6)
y_7 = clf.predict(X_test)

clf = clf.prune(2)
y_2 = clf.predict(X_test)

# Plot the results
import pylab as pl

pl.figure()
pl.scatter(X, y, c="k", label="data")
pl.plot(X_test, y_2, c="g", label="n_leaves=2", linewidth=2)
pl.plot(X_test, y_7, c="b", label="n_leaves=7", linewidth=2)
pl.plot(X_test, y_15, c="r", label="n_leaves=15", linewidth=2)
pl.xlabel("data")
pl.ylabel("target")
pl.title("Decision Tree Regression with levels of pruning")
pl.legend()
pl.show()
