"""
=====================================================
Gaussian process classification (GPC) on iris dataset
=====================================================

This example illustrates the predicted probability of GPC for an isotropic
and anisotropic RBF kernel on a two-dimensional version for the iris-dataset.
The anisotropic RBF kernel obtains slightly higher log-marginal-likelihood by
assigning different length-scales to the two feature dimensions.

"""

# Authors: The scikit-learn developers
# SPDX-License-Identifier: BSD-3-Clause

import matplotlib.pyplot as plt
import numpy as np

from sklearn import datasets
from sklearn.gaussian_process import GaussianProcessClassifier
from sklearn.gaussian_process.kernels import RBF
from sklearn.inspection import DecisionBoundaryDisplay

# import some data to play with
iris = datasets.load_iris()
X = iris.data[:, :2]  # we only take the first two features.
y = np.array(iris.target, dtype=int)

kernel = 1.0 * RBF([1.0])
gpc_rbf_isotropic = GaussianProcessClassifier(kernel=kernel).fit(X, y)
kernel = 1.0 * RBF([1.0, 1.0])
gpc_rbf_anisotropic = GaussianProcessClassifier(kernel=kernel).fit(X, y)

titles = ["Isotropic RBF", "Anisotropic RBF"]
fig, axes = plt.subplots(1, 2, figsize=(10, 5))
for ax, clf, title in zip(axes, (gpc_rbf_isotropic, gpc_rbf_anisotropic), titles):
    DecisionBoundaryDisplay.from_estimator(
        clf,
        X,
        response_method="predict",
        plot_method="contourf",
        grid_resolution=100,
        eps=0.5,
        ax=ax,
    )
    ax.scatter(X[:, 0], X[:, 1], c=np.array(["r", "g", "b"])[y], edgecolors=(0, 0, 0))
    ax.set_xlabel("Sepal length")
    ax.set_ylabel("Sepal width")
    ax.set_xticks(())
    ax.set_yticks(())
    ax.set_title(
        "%s, LML: %.3f" % (title, clf.log_marginal_likelihood(clf.kernel_.theta))
    )

plt.tight_layout()
plt.show()
