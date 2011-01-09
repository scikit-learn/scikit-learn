"""
===========================================================
A demo of feature_agglomeration on Lena
===========================================================

Author : Vincent Michel, 2010

"""
print __doc__

import time as time
import numpy as np
import scipy as sp
import pylab as pl
from scikits.learn.feature_extraction.image import img_to_graph
from scikits.learn.cluster import Ward, KMeans

###############################################################################
# Generate data
lena = sp.lena()
# Downsample the image by a factor of 4
lena = lena[::2, ::2] + lena[1::2, ::2] + lena[::2, 1::2] + lena[1::2, 1::2]
lena = lena[::2, ::2] + lena[1::2, ::2] + lena[::2, 1::2] + lena[1::2, 1::2]
mask = np.ones_like(lena).astype(bool)
X = np.atleast_2d(lena[mask])

###############################################################################
# Define the structure A of the data. Here a 10 nearest neighbors
adjacency_matrix = img_to_graph(mask, mask)

###############################################################################
# Compute feature agglomeration
k = 50
print "Compute feature agglomeration using Ward algorithm..."
ward = Ward(k=k).fit(X.T, adjacency_matrix=adjacency_matrix)
ward_label = np.reshape(ward.labels_, mask.shape)
ward_agglo = ward.inverse_transform(ward.transform(X.T))
print "Compute feature agglomeration using KMeans algorithm..."
kmeans = KMeans(k=k).fit(X.T)
kmeans_label = np.reshape(kmeans.labels_, mask.shape)
kmeans_agglo = kmeans.inverse_transform(kmeans.transform(X.T))



###############################################################################
# Plot the results on an image
pl.figure(figsize=(5, 5))
pl.imshow(np.reshape(ward_agglo, mask.shape),   cmap=pl.cm.gray)
for l in range(k):
    pl.contour(ward_label == l, contours=1,
            colors=[pl.cm.spectral(l/float(k)), ])
pl.xticks(())
pl.yticks(())
pl.figure(figsize=(5, 5))
pl.imshow(np.reshape(kmeans_agglo, mask.shape),   cmap=pl.cm.gray)
pl.xticks(())
pl.yticks(())
pl.show()
