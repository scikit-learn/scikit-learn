"""
==================================================
Demostrating `MDLPDiscretizer` in the iris dataset
==================================================

This example demonstrates the discretization performed
by the :class:`~sklearn.preprocessing.MDLPDiscretizer`
in the iris dataset.

This plot shows that the features where discretized
according with the (sub)optimal number of bins.
"""

# Author: Juan Carlos Alfaro Jiménez <JuanCarlos.Alfaro@uclm.es>

# License: BSD

from itertools import product

import numpy as np
import matplotlib.pyplot as plt
from sklearn.datasets import load_iris
from sklearn.experimental import enable_mdlp_discretizer  # noqa
from sklearn.preprocessing import MDLPDiscretizer

print(__doc__)

# Load the data
data = load_iris()
X, y = data.data, data.target
features = data.feature_names

# Get the output information
n_features = X.shape[1]
n_rows = n_features
n_cols = n_features

combinations = product(range(n_features), repeat=2)

# Discretizer for the iris dataset
discretizer = MDLPDiscretizer(encode="ordinal")

fig, axs = plt.subplots(nrows=n_rows,
                        ncols=n_cols,
                        figsize=(14, 9))

# Step size for the mesh
h = 0.02

# Blending value for the contour
alpha = 0.5

# Plot the discretization for each pair of features
for idx, (i, j) in enumerate(combinations):
    row = idx // n_cols
    col = idx % n_cols
    ax = axs[row, col]

    xi = X[:, i]
    xj = X[:, j]
    indices = [i, j]

    xi_min, xi_max = np.min(xi), np.max(xi)
    xj_min, xj_max = np.min(xj), np.max(xj)

    xx, yy = np.meshgrid(np.arange(xi_min, xi_max, h),
                         np.arange(xj_min, xj_max, h))

    grid = np.c_[np.ravel(xx), np.ravel(yy)]

    # Transform the dataset
    discretizer = discretizer.fit(X[:, indices], y)
    grid_encoded = discretizer.transform(grid)

    # Plot the horizontal and vertical stripes
    horizontal = grid_encoded[:, 0].reshape(xx.shape)
    vertical = grid_encoded[:, 1].reshape(yy.shape)

    ax.contourf(xx, yy, horizontal, alpha=alpha)
    ax.contourf(xx, yy, vertical, alpha=alpha)

    # Plot the data points
    ax.scatter(xi, xj, c=y, edgecolors="k")

    # Configure the axes for proper rendering
    ax.set_xlim(xi_min, xi_max)
    ax.set_ylim(xj_min, xj_max)

    if (row + 1) == n_rows:
        # Set the label for the x-axis
        ax.set_xlabel(features[col])
    if col == 0:
        # Set the label for the y-axis
        ax.set_ylabel(features[row])

# Adjust subplot to fit the figure area
plt.tight_layout()
plt.show()
