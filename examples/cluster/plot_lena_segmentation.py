"""
=========================================
Segmenting the picture of Lena in regions
=========================================

This example uses :ref:`spectral_clustering` on a graph created from
voxel-to-voxel difference on an image to break this image into multiple
partly-homogenous regions.

This procedure (spectral clustering on an image) is an efficient
approximate solution for finding normalized graph cuts.

Spectral clustering with kmeans will cluster samples in the embedding space 
with kmeans whereas discrete will iteratively search for the closest partition 
space to the embedding space.
"""
print __doc__

# Author: Gael Varoquaux <gael.varoquaux@normalesup.org>, Brian Cheung
# License: BSD

import time

import numpy as np
import scipy as sp
import pylab as pl

from sklearn.feature_extraction import image
from sklearn.cluster import spectral_clustering

lena = sp.misc.lena()
# Downsample the image by a factor of 4
lena = lena[::2, ::2] + lena[1::2, ::2] + lena[::2, 1::2] + lena[1::2, 1::2]
lena = lena[::2, ::2] + lena[1::2, ::2] + lena[::2, 1::2] + lena[1::2, 1::2]

# Convert the image into a graph with the value of the gradient on the
# edges.
graph = image.img_to_graph(lena)

# Take a decreasing function of the gradient: an exponential
# The smaller beta is, the more independent the segmentation is of the
# actual image. For beta=1, the segmentation is close to a voronoi
beta = 5
eps = 1e-6
graph.data = np.exp(-beta * graph.data / lena.std()) + eps

# Apply spectral clustering (this step goes much faster if you have pyamg
# installed)
N_REGIONS = 11
embedding_solvers = ['kmeans', 'discrete']
###############################################################################
# Visualize the resulting regions
pl.figure(figsize=(5, 10))
plot_num = 1
for embed_solve in embedding_solvers:
    t0 = time.time()
    labels = spectral_clustering(graph, n_clusters=N_REGIONS, embed_solve=embed_solve)
    t1 = time.time()
    labels = labels.reshape(lena.shape)

    pl.subplot(2,1,plot_num)
    pl.imshow(lena,   cmap=pl.cm.gray)
    for l in range(N_REGIONS):
        pl.contour(labels == l, contours=1,
                colors=[pl.cm.spectral(l / float(N_REGIONS)), ])
    pl.xticks(())
    pl.yticks(())
    pl.title('%s, %.2fs' % (embed_solve, (t1 - t0)))
        
    plot_num += 1
        
pl.show()