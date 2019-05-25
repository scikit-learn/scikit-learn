"""
================================================
Segmenting the picture of greek coins in regions
================================================

This example uses :ref:`spectral_clustering` on a graph created from
voxel-to-voxel difference on an image to break this image into multiple
partly-homogeneous regions.

This procedure (spectral clustering on an image) is an efficient
approximate solution for finding normalized graph cuts.

There are three options to assign labels:

* with 'kmeans' spectral clustering will cluster samples in the embedding space
  using a kmeans algorithm,
* with 'clusterQR' will cluster samples in the embedding space
  using a clusterQR algorithm,
* whereas 'discrete' will iteratively search for the closest partition
  space to the embedding space.

"""
print(__doc__)

# Author: Gael Varoquaux <gael.varoquaux@normalesup.org>, Brian Cheung
# Andrew Knyazev added clusterQR
# License: BSD 3 clause

import time

import numpy as np
from distutils.version import LooseVersion
from scipy.ndimage.filters import gaussian_filter
import matplotlib.pyplot as plt
import skimage
from skimage.data import coins
from skimage.transform import rescale

from sklearn.feature_extraction import image
from sklearn.cluster import spectral_clustering

# these were introduced in skimage-0.14
if LooseVersion(skimage.__version__) >= '0.14':
    rescale_params = {'anti_aliasing': False, 'multichannel': False}
else:
    rescale_params = {}

# load the coins as a numpy array
orig_coins = coins()

# Resize it to 20% of the original size to speed up the processing
# Applying a Gaussian filter for smoothing prior to down-scaling
# reduces aliasing artifacts.
smoothened_coins = gaussian_filter(orig_coins, sigma=2)
rescaled_coins = rescale(smoothened_coins, 0.2, mode="reflect",
                         **rescale_params)

# Convert the image into a graph with the value of the gradient on the
# edges.
graph = image.img_to_graph(rescaled_coins)

# Take a decreasing function of the gradient: an exponential
# The smaller beta is, the more independent the segmentation is of the
# actual image. For beta=1, the segmentation is close to a voronoi
beta = 10
eps = 1e-6
graph.data = np.exp(-beta * graph.data / graph.data.std()) + eps

# The actual number of regions in this example is 27: background and 26 coins
N_REGIONS = 26

#############################################################################
# Compute and visualize the resulting regions

# Any eigen_solver: 'arpack', 'lobpcg', 'amg' can be used. AMG is usually best
# It often helps the spectral clustering to compute a few extra eigenvectors
N_REGIONS_PLUS = 3

for assign_labels in ('kmeans', 'discretize', 'clusterQR'):
    t0 = time.time()
    labels = spectral_clustering(graph,
                                 n_clusters=(N_REGIONS + N_REGIONS_PLUS),
                                 assign_labels=assign_labels, random_state=42,
                                 eigen_solver='arpack')
    t1 = time.time()
    labels = labels.reshape(rescaled_coins.shape)

    plt.figure(figsize=(5, 5))
    plt.imshow(rescaled_coins, cmap=plt.get_cmap('gray'))
    plt.xticks(())
    plt.yticks(())
    title = 'Spectral clustering: %s, %.2fs' % (assign_labels, (t1 - t0))
    print(title)
    plt.title(title)
    for l in range(N_REGIONS):
        plt.contour(labels == l,
                    colors=[plt.cm.nipy_spectral((l+3) / float(N_REGIONS+3))])
        plt.pause(0.5)
plt.show()
