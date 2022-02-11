"""
============================
Faces dataset decompositions
============================

This example applies to :ref:`olivetti_faces_dataset` different unsupervised
matrix decomposition (dimension reduction) methods from the module
:py:mod:`sklearn.decomposition` (see the documentation chapter
:ref:`decompositions`) .


- Authors: Vlad Niculae, Alexandre Gramfort
- License: BSD 3 clause
"""

# %%
# Dataset preparation
# ---------------------------------------------------
#
# Loading and preprocessing the Olivetti faces dataset.

import logging
from time import time

from numpy.random import RandomState
import matplotlib.pyplot as plt

from sklearn.datasets import fetch_olivetti_faces
from sklearn.cluster import MiniBatchKMeans
from sklearn import decomposition

rng = RandomState(0)

# Display progress logs on stdout
logging.basicConfig(level=logging.INFO, format="%(asctime)s %(levelname)s %(message)s")

faces, _ = fetch_olivetti_faces(return_X_y=True, shuffle=True, random_state=rng)
n_samples, n_features = faces.shape

# Global centering (focus on one feature, centering all samples)
faces_centered = faces - faces.mean(axis=0)

# Local centering (focus on one sample, centering all features)
faces_centered -= faces_centered.mean(axis=1).reshape(n_samples, -1)

print("Dataset consists of %d faces" % n_samples)


# %%
# Define a base function to plot the gallery of faces.

n_row, n_col = 2, 3
n_components = n_row * n_col
image_shape = (64, 64)


def plot_gallery(title, images, n_col=n_col, n_row=n_row, cmap=plt.cm.gray):
    if n_col * n_row == 1:
        figsize = (3.5, 3.6)
    else:
        figsize = (2.0 * n_col, 2.26 * n_row)

    fig = plt.figure(figsize=figsize, facecolor="white", tight_layout={"pad": 0.3})
    fig.suptitle(title, size=16, wrap=True)
    for i, comp in enumerate(images):
        plt.subplot(n_row, n_col, i + 1)
        vmax = max(comp.max(), -comp.min())
        plt.imshow(
            comp.reshape(image_shape),
            cmap=cmap,
            interpolation="nearest",
            vmin=-vmax,
            vmax=vmax,
        )
        plt.xticks(())
        plt.yticks(())


# %%
# Let’s take a look at our data. Gray color indicates negative values,
# white indicates positive values.

plot_gallery("Faces from dataset", faces_centered[:n_components])

# %%
# Decomposition
# ---------------------------------------------------
#
# Initialise different estimators for decomposition and organise them all
# in a simple structure. Each tuple has the name and instance of the
# estimator and boolean flag if this estimator will be used
# on centred (preprocessed) data.

estimators = [
    (
        "Eigenfaces - PCA using randomized SVD",
        decomposition.PCA(
            n_components=n_components, svd_solver="randomized", whiten=True
        ),
        True,
    ),
    (
        "Non-negative components - NMF",
        decomposition.NMF(n_components=n_components, tol=5e-3),
        False,
    ),
    (
        "Independent components - FastICA",
        decomposition.FastICA(n_components=n_components, whiten=True),
        True,
    ),
    (
        "Sparse comp. - MiniBatchSparsePCA",
        decomposition.MiniBatchSparsePCA(
            n_components=n_components,
            alpha=0.8,
            n_iter=100,
            batch_size=3,
            random_state=rng,
        ),
        True,
    ),
    (
        "MiniBatchDictionaryLearning",
        decomposition.MiniBatchDictionaryLearning(
            n_components=15, alpha=0.1, n_iter=50, batch_size=3, random_state=rng
        ),
        True,
    ),
    (
        "Cluster centers - MiniBatchKMeans",
        MiniBatchKMeans(
            n_clusters=n_components,
            tol=1e-3,
            batch_size=20,
            max_iter=50,
            random_state=rng,
        ),
        True,
    ),
    (
        "Factor Analysis components - FA",
        decomposition.FactorAnalysis(n_components=n_components, max_iter=20),
        True,
    ),
]


# %%
# Next, fit each estimator on all images and plot some results.
# Each estimator extracts 6 components as vectors
# :math:`h \in \mathbf{R}^{4096}`. We just displayed these vectors
# in human-friendly visualisation as 64x64 pixel images.
#
# .. note::
#
#     The Eigenfaces estimator, via the :py:mod:`sklearn.decomposition.PCA`,
#     also provides a scalar `noise_variance_` (the mean of pixelwise variance)
#     that cannot be displayed as an image.

for name, estimator, center in estimators:
    print("- Extracting the top %d '%s'..." % (n_components, name))
    t0 = time()
    if center:
        data = faces_centered
    else:
        data = faces
    estimator.fit(data)
    train_time = time() - t0
    print("  done in %0.3fs" % train_time)
    if hasattr(estimator, "cluster_centers_"):
        components_ = estimator.cluster_centers_
    else:
        components_ = estimator.components_

    if (
        hasattr(estimator, "noise_variance_") and estimator.noise_variance_.ndim > 0
    ):  # skip the Eigenfaces case
        plot_gallery(
            "Pixelwise variance from \n %s" % name,
            estimator.noise_variance_.reshape(1, -1),
            n_col=1,
            n_row=1,
        )

    plot_gallery(name, components_[:n_components])

plt.show()
# %%
# Decomposition: Dictionary learning
# ---------------------------------------------------
#
# In the further section, let's consider :ref:`DictionaryLearning` more precisely.
# Dictionary learning is a problem that amounts to finding a sparse representation
# of the input data as a combination of simple elements. These simple elements form
# a dictionary. It is possible to constrain the dictionary and/or code to be positive
# to match constraints that may be present in the data.

dict_estimators = [
    (
        "Dictionary learning",
        decomposition.MiniBatchDictionaryLearning(
            n_components=15, alpha=0.1, n_iter=50, batch_size=3, random_state=rng
        ),
    ),
    (
        "Dictionary learning - positive dictionary",
        decomposition.MiniBatchDictionaryLearning(
            n_components=15,
            alpha=0.1,
            n_iter=50,
            batch_size=3,
            random_state=rng,
            positive_dict=True,
        ),
    ),
    (
        "Dictionary learning - positive code",
        decomposition.MiniBatchDictionaryLearning(
            n_components=15,
            alpha=0.1,
            n_iter=50,
            batch_size=3,
            fit_algorithm="cd",
            random_state=rng,
            positive_code=True,
        ),
    ),
    (
        "Dictionary learning - positive dictionary & code",
        decomposition.MiniBatchDictionaryLearning(
            n_components=15,
            alpha=0.1,
            n_iter=50,
            batch_size=3,
            fit_algorithm="cd",
            random_state=rng,
            positive_dict=True,
            positive_code=True,
        ),
    ),
]


# %%
# Plot the same samples from our dataset but with another colormap.
# Red indicates negative values, blue indicates positive values,
# and white represents zeros.

plot_gallery("Faces from dataset", faces_centered[:n_components], cmap=plt.cm.RdBu)

# %%
# Similar to the previous example, we train each estimator on all images.
#
# The results of the different positivity constraints are presented below.

for name, estimator in dict_estimators:
    print("- Extracting the top %d '%s'..." % (n_components, name))
    t0 = time()
    data = faces_centered
    estimator.fit(data)
    train_time = time() - t0
    print("  done in %0.3fs" % train_time)
    plot_gallery(name, estimator.components_[:n_components], cmap=plt.cm.RdBu)

plt.show()
