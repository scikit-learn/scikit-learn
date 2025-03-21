"""
=========================================================================
SVM: Effect of Regularization (C) on Maximum Margin Separating Hyperplane
========================================================================

This script demonstrates the concept of the maximum margin separating hyperplane in a two-class separable dataset using a Support Vector Machine (SVM) with a linear kernel and how different values of `C` influence the margin width.

- **Small C (e.g., 0.05)**: The model allows some misclassifications, resulting in a wider margin.
- **Moderate C (e.g., 1)**: The model balances classification accuracy and margin width.
- **Large C (e.g., 1000)**: The model prioritizes classifying all points correctly, leading to a narrower margin.

"""

# Authors: The scikit-learn developers
# SPDX-License-Identifier: BSD-3-Clause

import numpy as np
import matplotlib.pyplot as plt
from sklearn import svm
from sklearn.datasets import make_blobs
from sklearn.inspection import DecisionBoundaryDisplay

# create 40 separable points
X, y = make_blobs(n_samples=40, centers=2, cluster_std=1.5, random_state=6)

# define different values of C to observe its effect on the margin
C_values = [0.05, 1, 1000]

plt.figure(figsize=(15, 5))
for i, C_val in enumerate(C_values, 1):
    clf = svm.SVC(kernel="linear", C=C_val)
    clf.fit(X, y)

    plt.subplot(1, 3, i)
    plt.scatter(X[:, 0], X[:, 1], c=y, s=30, cmap=plt.cm.Paired, edgecolors="k")

    # plot the decision function
    ax = plt.gca()
    DecisionBoundaryDisplay.from_estimator(
        clf, X, plot_method="contour", colors="k", levels=[-1, 0, 1], alpha=0.5, linestyles=["--", "-", "--"], ax=ax
    )

    # plot support vectors
    ax.scatter(
        clf.support_vectors_[:, 0], clf.support_vectors_[:, 1],
        s=100, linewidth=1.5, facecolors="none", edgecolors="r", label="Support Vectors"
    )

    plt.title(f"SVM Decision Boundary (C={C_val})")
    plt.legend()

plt.show()
