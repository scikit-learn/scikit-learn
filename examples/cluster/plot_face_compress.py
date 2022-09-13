# -*- coding: utf-8 -*-
"""
=========================================================
Vector Quantization Example
=========================================================

Face, a 1024 x 768 size image of a raccoon face,
is used here to illustrate how KBinsDiscretizer is
used for vector quantization.

"""

# Code source: Gaël Varoquaux
# Modified for documentation by Jaques Grobler
# License: BSD 3 clause

import numpy as np
import scipy as sp
import matplotlib.pyplot as plt

from sklearn import preprocessing


try:  # SciPy >= 0.16 have face in misc
    from scipy.misc import face

    face = face(gray=True)
except ImportError:
    face = sp.face(gray=True)

n_bins = 5
np.random.seed(0)

X = face.reshape((-1, 1))  # We need an (n_sample, n_feature) array
est = preprocessing.KBinsDiscretizer(
    n_bins=n_bins, strategy="kmeans", encode="ordinal", random_state=0
)
est.fit(X)
values = est.bin_edges_[0]
labels = est.transform(X)

face_compressed = est.inverse_transform(labels)
face_compressed.shape = face.shape

vmin = face.min()
vmax = face.max()

# original face
plt.figure(1, figsize=(3, 2.2))
plt.title("original image")
plt.imshow(face, cmap=plt.cm.gray, vmin=vmin, vmax=256)

# compressed face
plt.figure(2, figsize=(3, 2.2))
plt.title("quantized image (kmeans)")
plt.imshow(face_compressed, cmap=plt.cm.gray, vmin=vmin, vmax=vmax)

# equal bins face
regular_values = np.linspace(0, 256, n_bins + 1)
regular_labels = np.searchsorted(regular_values, face) - 1
regular_values = 0.5 * (regular_values[1:] + regular_values[:-1])  # mean
regular_face = np.choose(regular_labels.ravel(), regular_values, mode="clip")
regular_face.shape = face.shape
plt.figure(3, figsize=(3, 2.2))
plt.title("quantized image (equal bins)")
plt.imshow(regular_face, cmap=plt.cm.gray, vmin=vmin, vmax=vmax)

# histogram
plt.figure(4, figsize=(3, 2.2))
plt.clf()
plt.axes([0.01, 0.01, 0.98, 0.98])
plt.hist(X, bins=256, color=".5", edgecolor=".5")
plt.yticks(())
plt.xticks(regular_values)
plt.title("image intensity histogram")
values = np.sort(values)
for center_1, center_2 in zip(values[:-1], values[1:]):
    plt.axvline(0.5 * (center_1 + center_2), color="b")

for center_1, center_2 in zip(regular_values[:-1], regular_values[1:]):
    plt.axvline(0.5 * (center_1 + center_2), color="b", linestyle="--")

plt.show()
# The number of unique values in the quantized image is reducd to only 5 values.
print(f"face: number of unique values: {len(np.unique(face))}")
print(f"labels: number of unique values: {len(np.unique(labels))}")
# A note on in-memory image compression
# Let's have a look at the memory size of these images.
print(f"face bytes: {face.nbytes}")
print(f"labels bytes: {labels.nbytes}")
# Although the quantized image has far fewer unique values,
# its memory size is 4 times the original image.
# The reason behind this is the data types used.
print(f"face dtype: {face.dtype}")
print(f"labels dtype: {labels.dtype}")
# Let us cast the compressed image to the same type as the input
labels_uint8 = labels.astype(np.uint8)
# Now, both images have the same size
print(f"face bytes: {face.nbytes}")
print(f"labels bytes: {labels_uint8.nbytes}")
# The cluster labels are in the range 0,1,2,3,4
# While this could be compressed to a 4 bit integer,
# uint8 is as small as we can go with numpy.
