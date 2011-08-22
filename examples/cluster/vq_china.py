# -*- coding: utf-8 -*-
"""
=========================================
Vector Quantization of Lena using k-means
=========================================

Performs a Vector Quatization of an image, reducing the 
number of colors required to show the image.
"""
print __doc__

import numpy as np
import pylab as pl
from scikits.learn.cluster import KMeans
from scikits.learn.datasets import load_sample_images
# Try to import Image and imresize from PIL. We do this here to prevent
# this module from depending on PIL.
try:
    try:
        from scipy.misc import Image
    except ImportError:
        from scipy.misc.pilutil import Image
except ImportError:
    raise ImportError("The Python Imaging Library (PIL)"
                      "is required to load data from jpeg files")

# Get all sample images and obtain just china.jpg
sample_images = load_sample_images()
index = sample_images.filenames.index("china.jpg")
image_data = sample_images.images[index]

# Load Image and transform to a 2D numpy array.
w, h, d = original_shape = tuple(image_data.shape)
image_array = np.reshape(image_data, (w * h, 3))

# Take a sample of the data.
sample_indices = range(len(image_array))
np.random.shuffle(sample_indices)
sample_indices = sample_indices[:int(len(image_array) * 0.5)]
sample_data = image_array[sample_indices]

# Perform Vector Quantisation with 256 clusters.
k = 256
kmeans = KMeans(k=k)
kmeans.fit(image_array)
# Get labels for all points
labels = kmeans.predict(image_array)
# Save the reduced dataset. Only the centroids and labels need to be saved.
reduced_image = (kmeans.cluster_centers_, labels)

def recreate_image(centroids, labels, w, h):
    # Recreates the (compressed) image from centroids, labels and dimensions
    d = len(centroids[0])
    image = np.zeros((w, h, d))
    label_num = 0
    for i in range(w):
        for j in range(h):
            image[i][j] = centroids[labels[label_num]]
            print labels[label_num], label_num
            label_num += 1
    return image

# Display all results, alongside original image
pl.figure()
ax = pl.axes([0,0,1,1], frameon=False)
ax.set_axis_off()
centroids, labels = reduced_image
im = pl.imshow(recreate_image(centroids, labels, w, h))

show()

