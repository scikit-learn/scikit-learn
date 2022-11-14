"""
=======================
IsolationForest example
=======================

An example using :class:`~sklearn.ensemble.IsolationForest` for anomaly
detection.

The IsolationForest 'isolates' observations by randomly selecting a feature and
then randomly selecting a split value between the maximum and minimum values of
the selected feature.

Since recursive partitioning can be represented by a tree structure, the number
of splittings required to isolate a sample is equivalent to the path length from
the root node to the terminating node.

This path length, averaged over a forest of such random trees, is a measure of
normality and our decision function.

Random partitioning produces noticeable shorter paths for anomalies. Hence, when
a forest of random trees collectively produce shorter path lengths for
particular samples, they are highly likely to be anomalies.

"""

# %%
# Data generation
# ---------------
#
# We generate two clusters (each one containing `n_samples`) by randomly
# sampling the standard normal distribution as returned by `numpy.random.randn`.
# One of them is spherical and the other one is slightly deformed.
#
# For consistency with the :class:`~sklearn.ensemble.IsolationForest` notation,
# the inliers (i.e. the gaussian clusters) are assigned a ground truth label `1`
# whereas the outliers (created with `numpy.random.uniform`) are assigned the
# label `-1`.

import numpy as np
from sklearn.model_selection import train_test_split

n_samples = 120
n_outliers = 40
rng = np.random.RandomState(0)
C = np.array([[0.5, -0.1], [0.7, 0.4]])
cluster_1 = 0.4 * np.dot(rng.randn(n_samples, 2), C) + np.array([2, 2])  # general
cluster_2 = 0.3 * rng.randn(n_samples, 2) + np.array([-2, -2])  # spherical
outliers = rng.uniform(low=-4, high=4, size=(n_outliers, 2))

X = np.concatenate([cluster_1, cluster_2, outliers])
y = np.concatenate(
    [np.ones((2 * n_samples), dtype=int), -np.ones((n_outliers), dtype=int)]
)

X_train, X_test, y_train, y_test = train_test_split(X, y, stratify=y, random_state=42)

# %%
# We can visualize the resulting clusters:

import matplotlib.pyplot as plt

plt.scatter(X[:, 0], X[:, 1], c=y, s=20, edgecolor="k")
plt.title("Gaussian Mixture clusters")
plt.show()

# %%
# We train the model and use the class
# :class:`~sklearn.inspection.DecisionBoundaryDisplay` to visualize a discrete
# decision boundary to determine whether a particular sample is an outlier or
# not.

import matplotlib.pyplot as plt
from sklearn.inspection import DecisionBoundaryDisplay
from sklearn.ensemble import IsolationForest

clf = IsolationForest(max_samples=100, random_state=rng)
clf.fit(X_train)

disp = DecisionBoundaryDisplay.from_estimator(
    clf,
    X,
    response_method="predict",
    alpha=0.5,
)
scatter = disp.ax_.scatter(X[:, 0], X[:, 1], c=y, s=20, edgecolor="k")
disp.ax_.legend(
    *scatter.legend_elements(),
    title="True class",
    loc="upper left",
)
disp.ax_.set_title("Binary decision boundary of IsolationForest")
plt.show()

# %%
# By setting the `response_method="decision_function"`, the
# :class:`~sklearn.inspection.DecisionBoundaryDisplay` plots instead the measure
# of normality of an observation, which is given by the depth of the leaf (or
# equivalently the number of splits) required to isolate a sample in a given
# position.

disp = DecisionBoundaryDisplay.from_estimator(
    clf,
    X,
    response_method="decision_function",
    alpha=0.5,
)
scatter = disp.ax_.scatter(X[:, 0], X[:, 1], c=y, s=20, edgecolor="k")
disp.ax_.legend(
    *scatter.legend_elements(),
    title="True class",
    loc="upper left",
)
disp.ax_.set_title("Path length decision boundary of IsolationForest")
plt.show()
