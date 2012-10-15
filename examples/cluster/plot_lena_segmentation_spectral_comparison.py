"""
===============================================================
Spectral clustering embedding comparison for image segmentation
===============================================================

This example compares different methods to handle the embedding space of
:ref:`spectral_clustering`.

In these settings, the spectral clustering approach solves the problem
know as 'normalized graph cuts': the image is seen as a graph of
connected voxels, and the spectral clustering algorithm amounts to
choosing graph cuts defining regions while minimizing the ratio of the
gradient along the cut, and the volume of the region.

As the algorithm tries to balance the volume (ie balance the region
sizes), if we take circles with different sizes, the segmentation fails.

In addition, as there is no useful information in the intensity of the image,
or its gradient, we choose to perform the spectral clustering on a graph
that is only weakly informed by the gradient. This is close to performing
a Voronoi partition of the graph.

In addition, we use the mask of the objects to restrict the graph to the
outline of the objects. In this example, we are interested in
separating the objects one from the other, and not from the background.

Spectral clustering with kmeans will cluster samples in the embedding space 
with kmeans whereas discrete will iteratively search for the closest partition 
space to the embedding space.  This examples shows the properties of the
segmented regions for various number of clusters.  The computation time is
also shown to illustrate the scalability of the embedding solvers for the
specified number of clusters.
"""
print __doc__

# Author: Brian Cheung
# Modification of original lena segmentation example by
# Gael Varoquaux <gael.varoquaux@normalesup.org>
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

pl.figure(figsize=(5,15))

segmentations = [4,6,8,10,12]
embedding_solvers = ['kmeans', 'discrete']
plot_num = 1
for n_regions in segmentations:
    for embed_solve in embedding_solvers:
        t0 = time.time()
        labels = spectral_clustering(graph, n_clusters=n_regions, embed_solve=embed_solve)
        t1 = time.time()
        labels = labels.reshape(lena.shape)
        
        pl.subplot(len(segmentations), len(embedding_solvers), plot_num)
        
        pl.imshow(lena,   cmap=pl.cm.gray)
        for l in range(n_regions):
            pl.contour(labels == l, contours=1,
                    colors=[pl.cm.spectral(l / float(n_regions)), ])
        pl.xticks(())
        pl.yticks(())
        if plot_num <= len(embedding_solvers):
            pl.title('%s\nn_regions: %d, %.2fs' % (embed_solve, n_regions, (t1 - t0)))
        else:
            pl.title('n_regions: %d, %.2fs' % (n_regions, (t1 - t0)))
        plot_num += 1
    
pl.show()
