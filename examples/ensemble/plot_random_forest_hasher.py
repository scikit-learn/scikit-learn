"""
===================================================
Hashing feature transformation using Random Forests
===================================================

RandomForestHasher provide a way to map data to a
very high-dimensional, sparse representation, which might
be beneficial for classification.
The mapping is completely unsupervised and very efficient.

This example visualizes the partitionings given by several
trees and shows how the transformation can also be used for
non-linear dimensionality reduction or manifold learning.

Points that are neighboring often share the same leaf of a tree will share
large parts of their hashed representation. This allows to
separate two circles simply based on the principal components of the
transformed data.

In the high-dimensional space, a simple classifier if often
enough for a good fit. For sparse binary data, BernoulliNB
is particularly well-suited. The bottom row compares the
decision boundary obtained by BernoulliNB in the transformed
space with an ExtraTreesClassifier forests learned on the
original data.
"""
import pylab as pl
import numpy as np

from sklearn.datasets import make_circles
from sklearn.ensemble import RandomForestHasher, ExtraTreesClassifier
from sklearn.decomposition import RandomizedPCA
from sklearn.naive_bayes import BernoulliNB

# make a synthetic dataset
X, y = make_circles(factor=0.5, random_state=0, noise=0.05)

# use RandomForestHasher to transform data
hasher = RandomForestHasher(n_estimators=10, random_state=0, max_depth=3)
X_transformed = hasher.fit_transform(X)

# Visualize result using PCA
pca = RandomizedPCA(n_components=2)
X_reduced = pca.fit_transform(X_transformed)

# Learn a Naive Bayes classifier on the transformed data
nb = BernoulliNB()
nb.fit(X_transformed, y)


# Learn an ExtraTreesClassifier for comparison
trees = ExtraTreesClassifier(max_depth=3, n_estimators=10, random_state=0)
trees.fit(X, y)


# scatter plot of original and reduced data
fig, axes = pl.subplots(2, 2, figsize=(8, 8))
axes = axes.ravel()
axes[0].scatter(X[:, 0], X[:, 1], c=y, s=50)
axes[0].set_title("Original Data")
axes[1].scatter(X_reduced[:, 0], X_reduced[:, 1], c=y, s=50)
axes[1].set_title("PCA reduction of transformed data")

# Plot the decision in original space. For that, we will asign a color to each
# point in the mesh [x_min, m_max] x [y_min, y_max].
h = .01
x_min, x_max = X[:, 0].min() - .5, X[:, 0].max() + .5
y_min, y_max = X[:, 1].min() - .5, X[:, 1].max() + .5
xx, yy = np.meshgrid(np.arange(x_min, x_max, h),
        np.arange(y_min, y_max, h))

# transform grid using RandomForestHasher
transformed_grid = hasher.transform(np.c_[xx.ravel(), yy.ravel()])
y_grid_pred = nb.predict_proba(transformed_grid)[:, 1]
axes[2].set_title("Naive Bayes on Transformed data")
axes[2].pcolormesh(xx, yy, y_grid_pred.reshape(xx.shape))
axes[2].scatter(X[:, 0], X[:, 1], c=y, s=50)
axes[2].set_ylim(-1.4, 1.4)
axes[2].set_xlim(-1.4, 1.4)

# transform grid using RandomForestHasher
y_grid_pred = trees.predict_proba(np.c_[xx.ravel(), yy.ravel()])[:, 1]
axes[3].set_title("ExtraTrees predictions")
axes[3].pcolormesh(xx, yy, y_grid_pred.reshape(xx.shape))
axes[3].scatter(X[:, 0], X[:, 1], c=y, s=50)
axes[3].set_ylim(-1.4, 1.4)
axes[3].set_xlim(-1.4, 1.4)
fig.subplots_adjust(left=0.02, right=0.98, top=0.94, bottom=0.02)

for ax in axes:
    ax.set_xticks(())
    ax.set_yticks(())
pl.show()
