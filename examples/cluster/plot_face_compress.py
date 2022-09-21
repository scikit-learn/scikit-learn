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
import scipy as sp
import matplotlib.pyplot as plt

from sklearn import cluster


try:  # SciPy >= 0.16 have face in misc
    from scipy.misc import face

    face = face(gray=True)
except ImportError:
    face = sp.misc.face(gray=True)

N_CLUSTERS = 5
np.random.seed(0)

X = face.reshape((-1, 1))  # We need an (n_sample, n_feature) array
k_means = cluster.KMeans(n_clusters=N_CLUSTERS, n_init=4)
k_means.fit(X)
values = k_means.cluster_centers_.squeeze()
labels = k_means.labels_

# create an array from labels and values
face_compressed = np.choose(labels, values)
face_compressed.shape = face.shape

vmin = face.min()
vmax = face.max()

# original face
plt.figure(1, figsize=(3, 2.2))
plt.imshow(face, cmap=plt.cm.gray, vmin=vmin, vmax=256)

# compressed face
plt.figure(2, figsize=(3, 2.2))
plt.imshow(face_compressed, cmap=plt.cm.gray, vmin=vmin, vmax=vmax)

# equal bins face
regular_values = np.linspace(0, 256, N_CLUSTERS + 1)
regular_labels = np.searchsorted(regular_values, face) - 1
regular_values = 0.5 * (regular_values[1:] + regular_values[:-1])  # mean
regular_face = np.choose(regular_labels.ravel(), regular_values, mode="clip")
regular_face.shape = face.shape
plt.figure(3, figsize=(3, 2.2))
plt.imshow(regular_face, cmap=plt.cm.gray, vmin=vmin, vmax=vmax)

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

#Word Of Caution:

#Theoretically, the above example should demonstrate compression.
#But there's a catch when we try to check the size of images;
#the compressed image actually ends up taking more memory
#space than the image itself. Even though the number of
#unique values in the image is greater compared to the
#compressed image.The catch comes from the fact that
#numpy arrays are contiguous and are required to have
#preallocated memory for each object. This is where
#the numpy speed comes from. But for example, as above,
#it has caused problems. The kmeans does reduce the
#dimensionality or number of unique values required to
#represent the image,but they are still stored in the same length nparray.

#Colored Image version With actual Compression:

face = sp.misc.face(gray=False)

N_CLUSTERS = 64  # clusters of colors
np.random.seed(0)

plt.imshow(face)

w, h, d = original_shape = tuple(
    face.shape
)  # a colored image contains 3 stacked matrices

face = np.array(face, dtype=np.float64) / 255  # small number easier to multiply
face_array = np.reshape(face, (w * h, d))  # reshape the image for Kmeans
kmeans = cluster.KMeans(n_clusters=N_CLUSTERS, random_state=0).fit(face_array)
labels = kmeans.predict(face_array)


plt.figure(1)
plt.clf()
plt.axis("on")
plt.title("Original image")
plt.imshow(face)

plt.figure(2)
plt.clf()
plt.axis("on")
plt.title("Quantised image")
img_re = kmeans.cluster_centers_[labels].reshape(w, h, -1)
plt.imshow(img_re)

plt.figure(3)
plt.clf()
plt.axis("on")
plt.title("Quantised image Compressed")
img_re_compressed = (img_re * 255).astype("uint8")
plt.imshow(img_re_compressed)

# number of elements multiplied by byte size of each element
img_re_size = img_re.itemsize * img_re.size
img_re_c_size = img_re_compressed.itemsize * img_re_compressed.size
# number of unique elements reduce in compressed image due to clustering
face_uniq = len(np.unique(face))
img_re_c_uniq = len(np.unique(img_re_compressed))

print(f"Size of reconstructed image: {img_re_size} bytes")
print(f"Size of compressed reconstructed image: {img_re_c_size} bytes")
print(f"number of unique values in original image: {face_uniq}")
print(f"number of unique values in reconstructed image: {img_re_c_uniq}")
