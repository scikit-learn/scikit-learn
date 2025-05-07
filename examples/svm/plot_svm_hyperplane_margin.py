"""
=========================================================================
SVM: Effect of Regularization (C) on Maximum Margin Separating Hyperplane
=========================================================================

This script demonstrates the concept of maximum margin separating hyperplane
in a two-class separable dataset using a Support Vector Machine (SVM)
with a linear kernel and how different values of `C` influence margin width.

- **Small C (e.g., 0.05)**:
    - Allows some misclassifications, resulting in wider margin.
- **Moderate C (e.g., 1)**:
    - Balances classification accuracy and margin width.
- **Large C (e.g., 1000)**:
    - Prioritizes classifying all points correctly, leading to narrower margin.

"""

# Authors: The scikit-learn developers
# SPDX-License-Identifier: BSD-3-Clause

# %%
import matplotlib.pyplot as plt

from sklearn import svm
from sklearn.datasets import make_blobs
from sklearn.inspection import DecisionBoundaryDisplay

# %%
# Create 40 separable points
X, y = make_blobs(n_samples=40, centers=2, cluster_std=1.5, random_state=6)

# %%
# Define different values of C to observe its effect on the margin
C_values = [0.05, 1, 1000]

# %%
# Visualize
plt.figure(figsize=(12, 4))
for i, C_val in enumerate(C_values, 1):
    clf = svm.SVC(kernel="linear", C=C_val)
    clf.fit(X, y)
    y_pred = clf.predict(X)
    misclassified = y_pred != y

    plt.subplot(1, 3, i)
    plt.scatter(
        X[:, 0], X[:, 1],
        c=y,
        s=30,
        cmap=plt.cm.Paired,
        edgecolors="k"
    )
    # misclassified samples
    plt.scatter(
        X[misclassified, 0],
        X[misclassified, 1],
        facecolors="none",
        edgecolors="k",
        s=80,
        linewidths=1.5,
        label="Misclassified",
    )

    # plot the decision function
    ax = plt.gca()
    DecisionBoundaryDisplay.from_estimator(
        clf,
        X,
        plot_method="contour",
        colors="k",
        levels=[-1, 0, 1],
        alpha=0.5,
        linestyles=["--", "-", "--"],
        ax=ax,
    )

    # plot support vectors
    ax.scatter(
        clf.support_vectors_[:, 0],
        clf.support_vectors_[:, 1],
        s=120,
        linewidth=1.5,
        facecolors="none",
        edgecolors="r",
        label="Support Vectors",
    )

    plt.title(f"SVM Decision Boundary (C={C_val})")
    plt.xlabel("Feature 1")
    plt.ylabel("Feature 2")
    plt.legend()

plt.tight_layout()
plt.show()

# %% [markdown]
# - **Small `C` (e.g., 0.01, 0.05)**:
#   - Use when:
#     - You expect noisy or overlapping data.
#     - You can tolerate some misclassification in training.
#     - Your priority is better generalization on unseen data.
#   - Note:
#     - May underfit if the margin is too lenient.
# - **Moderate `C` (e.g., 1)**:
#   - Use when:
#     - You're unsure about noise levels.
#     - You want good balance between margin width and classification accuracy.
# - **Large `C` (e.g., 1000)**:
#   - Use when:
#     - The data is clean and linearly separable.
#     - You want to avoid any training misclassification.
#   - Note:
#     - May overfit noisy data by trying to classify all samples correctly.
