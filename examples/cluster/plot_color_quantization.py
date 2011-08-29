# -*- coding: utf-8 -*-
"""
============================================
Vector Quantization of a photo using k-means
============================================

Performs a pixel-wise Vector Quantization (VQ) of an image of the summer palace
(China), reducing the number of colors required to show the image from 96,615
unique colors to 256, while preserving the overall appearance quality.

In this example, pixels are represented in a 3D-space and K-means is used to
find 256 clusters. In the image processing literature, the codebook obtained
from K-means (the cluster centers) is called the color palette. Using a
256-color palette, pixels can be encoded with a single byte (the index of the
closest color in the palette) whereas an RGB encoding requires 3 bytes per
pixel. The GIF file format, for example, uses such a palette.
"""
print __doc__
import numpy as np
import pylab as pl
from scikits.learn.cluster import KMeans
from scikits.learn.datasets import load_sample_image
from scikits.learn.utils import shuffle
from time import time

# Load the Summer Palace photo
china = load_sample_image("china.jpg")

# Convert to floats instead of the default 8 bits integer coding. Dividing by
# 255 is important so that pl.imshow behaves works well on foat data (need to
# be in the range [0-1]
china = np.array(china, dtype=np.float64) / 255

# Load Image and transform to a 2D numpy array.
w, h, d = original_shape = tuple(china.shape)
assert d == 3
image_array = np.reshape(china, (w * h, d))

print "Fitting estimator on a small sub-sample of the data"
t0 = time()
image_array_sample = shuffle(image_array, random_state=0)[:1000]
kmeans = KMeans(k=256, random_state=0).fit(image_array_sample)
print "done in %0.3fs." % (time() - t0)

# Get labels for all points
print "Predicting labels on the full image"
t0 = time()
labels = kmeans.predict(image_array)
print "done in %0.3fs." % (time() - t0)


def recreate_image(codebook, labels, w, h):
    """Recreate the (compressed) image from the code book & labels"""
    d = codebook.shape[1]
    image = np.zeros((w, h, d))
    label_idx = 0
    for i in range(w):
        for j in range(h):
            image[i][j] = codebook[labels[label_idx]]
            label_idx += 1
    return image

# Display all results, alongside original image
pl.figure(1)
pl.clf()
ax = pl.axes([0, 0, 1, 1])
pl.axis('off')
pl.title('Original image (96,615 colors)')
pl.imshow(china)

pl.figure(2)
pl.clf()
ax = pl.axes([0, 0, 1, 1])
pl.axis('off')
pl.title('Quantized image (256 colors)')
pl.imshow(recreate_image(kmeans.cluster_centers_, labels, w, h))
pl.show()
