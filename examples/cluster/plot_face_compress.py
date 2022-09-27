# -*- coding: utf-8 -*-
"""
=========================================================
Vector Quantization Example
=========================================================

Face, a 1024 x 768 size image of a raccoon face,
is used here to illustrate how `k`-means is
used for vector quantization.

"""

# Code source: Gaël Varoquaux
# Modified for documentation by Jaques Grobler
# License: BSD 3 clause

import numpy as np
import matplotlib.pyplot as plt

from sklearn import cluster

try:  # Scipy >= 1.10
    from scipy.datasets import face
except ImportError:
    from scipy.misc import face

racoon_face = face(gray=True)

n_clusters = 5
np.random.seed(0)

X = racoon_face.reshape((-1, 1))  # We need an (n_sample, n_feature) array
k_means = cluster.KMeans(n_clusters=n_clusters, n_init=4)
k_means.fit(X)
values = k_means.cluster_centers_.squeeze()
labels = k_means.labels_

# create an array from labels and values
racoon_face_compressed = np.choose(labels, values)
racoon_face_compressed.shape = racoon_face.shape

vmin = racoon_face.min()
vmax = racoon_face.max()

# original racoon_face
plt.figure(1, figsize=(3, 2.2))
plt.imshow(racoon_face, cmap=plt.cm.gray, vmin=vmin, vmax=256)

# compressed racoon_face
plt.figure(2, figsize=(3, 2.2))
plt.imshow(racoon_face_compressed, cmap=plt.cm.gray, vmin=vmin, vmax=vmax)

# equal bins racoon_face
regular_values = np.linspace(0, 256, n_clusters + 1)
regular_labels = np.searchsorted(regular_values, racoon_face) - 1
regular_values = 0.5 * (regular_values[1:] + regular_values[:-1])  # mean
regular_racoon_face = np.choose(regular_labels.ravel(), regular_values, mode="clip")
regular_racoon_face.shape = racoon_face.shape
plt.figure(3, figsize=(3, 2.2))
plt.imshow(regular_racoon_face, cmap=plt.cm.gray, vmin=vmin, vmax=vmax)

# histogram
plt.figure(4, figsize=(3, 2.2))
plt.clf()
plt.axes([0.01, 0.01, 0.98, 0.98])
plt.hist(X, bins=256, color=".5", edgecolor=".5")
plt.yticks(())
plt.xticks(regular_values)
values = np.sort(values)
for center_1, center_2 in zip(values[:-1], values[1:]):
    plt.axvline(0.5 * (center_1 + center_2), color="b")

for center_1, center_2 in zip(regular_values[:-1], regular_values[1:]):
    plt.axvline(0.5 * (center_1 + center_2), color="b", linestyle="--")

plt.show()
