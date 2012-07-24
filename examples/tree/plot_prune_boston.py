"""
============================================
Cross validated scores of the boston dataset
============================================

"""
print __doc__

import numpy as np
from sklearn.datasets import load_boston
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


boston = load_boston()
clf = tree.DecisionTreeRegressor(max_depth=8)

#Compute the cross validated scores
scores = tree.prune_path(clf, boston.data, boston.target,
                                    max_n_leaves=20, n_iterations=10,
                                    random_state=0)

plot_pruned_path(scores)
