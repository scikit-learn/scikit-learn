"""
===================================================
Hashing feature transformation using Random Forests
===================================================

RandomForestEmbedding provide a way to map data to a
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
from sklearn.ensemble import RandomForestEmbedding, ExtraTreesClassifier
from sklearn.decomposition import RandomizedPCA
from sklearn.naive_bayes import BernoulliNB

# make a synthetic dataset
X, y = make_circles(factor=0.5, random_state=0, noise=0.05)

# use RandomForestEmbedding to transform data
hasher = RandomForestEmbedding(n_estimators=10, random_state=0, max_depth=3)
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
fig = pl.figure(figsize=(9, 8))

ax = pl.subplot(221)
ax.scatter(X[:, 0], X[:, 1], c=y, s=50)
ax.set_title("Original Data")
ax.set_xticks(())
ax.set_yticks(())

ax = pl.subplot(222)
ax.scatter(X_reduced[:, 0], X_reduced[:, 1], c=y, s=50)
ax.set_title("PCA reduction of transformed data")
ax.set_xticks(())
ax.set_yticks(())

# Plot the decision in original space. For that, we will asign a color to each
# point in the mesh [x_min, m_max] x [y_min, y_max].
h = .01
x_min, x_max = X[:, 0].min() - .5, X[:, 0].max() + .5
y_min, y_max = X[:, 1].min() - .5, X[:, 1].max() + .5
xx, yy = np.meshgrid(np.arange(x_min, x_max, h),
        np.arange(y_min, y_max, h))

# transform grid using RandomForestEmbedding
transformed_grid = hasher.transform(np.c_[xx.ravel(), yy.ravel()])
y_grid_pred = nb.predict_proba(transformed_grid)[:, 1]

ax = pl.subplot(223)
ax.set_title("Naive Bayes on Transformed data")
ax.pcolormesh(xx, yy, y_grid_pred.reshape(xx.shape))
ax.scatter(X[:, 0], X[:, 1], c=y, s=50)
ax.set_ylim(-1.4, 1.4)
ax.set_xlim(-1.4, 1.4)
ax.set_xticks(())
ax.set_yticks(())

# transform grid using ExtraTreesClassifier
y_grid_pred = trees.predict_proba(np.c_[xx.ravel(), yy.ravel()])[:, 1]

ax = pl.subplot(224)
ax.set_title("ExtraTrees predictions")
ax.pcolormesh(xx, yy, y_grid_pred.reshape(xx.shape))
ax.scatter(X[:, 0], X[:, 1], c=y, s=50)
ax.set_ylim(-1.4, 1.4)
ax.set_xlim(-1.4, 1.4)
ax.set_xticks(())
ax.set_yticks(())

pl.tight_layout()
pl.show()
