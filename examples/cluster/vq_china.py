# -*- coding: utf-8 -*-
"""
=========================================
Vector Quantization of Lena using k-means
=========================================

Performs a Vector Quantization of an image, reducing the number of colors
required to show the image.
"""
print __doc__
import numpy as np
import pylab as pl
from scikits.learn.preprocessing import Scaler
from scikits.learn.cluster import KMeans
from scikits.learn.datasets import load_sample_image
from scikits.learn.utils import shuffle
from time import time

# Get all sample images and obtain just china.jpg
china = load_sample_image("china.jpg")

# Convert to floats instead of the default 8 bits integer coding
china = np.array(china, dtype=np.float64)

# Load Image and transform to a 2D numpy array.
w, h, d = original_shape = tuple(china.shape)
assert d == 3
image_array = np.reshape(china, (w * h, d))
scaler = Scaler()
image_scaled = scaler.fit_transform(image_array)

print "Fitting estimator on a small sub-sample of the data"
t0 = time()
image_scaled_sample = shuffle(image_scaled, random_state=0)[:1000]
kmeans = KMeans(k=10, max_iter=1000).fit(image_scaled_sample)
print "done in %0.3fs." % (time() - t0)

# Get labels for all points
print "Predicting labels on the full image"
t0 = time()
labels = kmeans.predict(image_scaled)
print "done in %0.3fs." % (time() - t0)

def recreate_image(codebook, labels, w, h):
    # Recreates the (compressed) image from the code book, labels and dimensions
    d = codebook.shape[1]
    image = np.zeros((w, h, d))
    label_idx = 0
    for i in range(w):
        for j in range(h):
            image[i][j] = codebook[labels[label_idx]]
            label_idx += 1
    return image

# Display all results, alongside original image
pl.figure()
ax = pl.axes([0, 0, 1, 1], frameon=False)
ax.set_axis_off()
pl.imshow(china)

pl.figure()
ax = pl.axes([0, 0, 1, 1], frameon=False)
ax.set_axis_off()

# transform back the images in the original space
codebook = scaler.inverse_transform(kmeans.cluster_centers_)

pl.imshow(recreate_image(kmeans.cluster_centers_, labels, w, h))
pl.show()
