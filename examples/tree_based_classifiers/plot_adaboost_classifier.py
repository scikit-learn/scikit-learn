"""
=========================================================
AdaBoost on Gaussian Quantiles
=========================================================

This example fits an :class:`~sklearn.ensemble.AdaBoostClassifier` on a
dataset consisting of Gaussian Quantiles, which is non-linearly separable.
"""

# Authors: The scikit-learn developers
# SPDX-License-Identifier: BSD-3-Clause

import matplotlib.pyplot as plt
import numpy as np
from sklearn.datasets import make_gaussian_quantiles
from sklearn.ensemble import AdaBoostClassifier
from sklearn.tree import DecisionTreeClassifier
from sklearn.inspection import DecisionBoundaryDisplay

# %%
# Generate Data
# -------------
X1, y1 = make_gaussian_quantiles(
    cov=2.0, n_samples=200, n_features=2, n_classes=2, random_state=1
)
X2, y2 = make_gaussian_quantiles(
    mean=(3, 3), cov=1.5, n_samples=300, n_features=2, n_classes=2, random_state=1
)
X = np.concatenate((X1, X2))
y = np.concatenate((y1, -y2 + 1))

# %%
# Train AdaBoost
# --------------
# We use a DecisionTreeClassifier with max_depth=1 (Decision Stump) as the base estimator.
# We compare two AdaBoost classifiers with different numbers of estimators.
bdt = AdaBoostClassifier(
    DecisionTreeClassifier(max_depth=1), algorithm="SAMME", n_estimators=200
)
bdt.fit(X, y)

# %%
# Plot Decision Boundary
# ----------------------
plot_colors = "br"
plot_step = 0.02
class_names = "AB"

plt.figure(figsize=(10, 5))

# Plot the decision boundaries
ax = plt.subplot(111)
DecisionBoundaryDisplay.from_estimator(
    bdt,
    X,
    cmap=plt.cm.Paired,
    response_method="predict",
    ax=ax,
    shading="auto",
)

# Plot the training points
for i, n, c in zip(range(2), class_names, plot_colors):
    idx = np.where(y == i)
    plt.scatter(
        X[idx, 0],
        X[idx, 1],
        c=c,
        cmap=plt.cm.Paired,
        s=20,
        edgecolor="k",
        label=f"Class {n}",
    )
plt.xlim([-5, 5])
plt.ylim([-5, 5])
plt.legend(loc="upper right")
plt.xlabel("x")
plt.ylabel("y")
plt.title("Decision Boundary of AdaBoostClassifier")

plt.show()
