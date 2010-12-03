"""
===========================================================
A demo of structure ward on Lena
===========================================================

Author : Vincent Michel, 2010

"""
print __doc__


import time as time
import numpy as np
import scipy as sp
import pylab as pl
from scikits.learn.feature_extraction.image import img_to_graph
from scikits.learn.cluster import Ward, plot_dendrogram

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
adjacency_matrix = img_to_graph(mask,mask)

###############################################################################
# Compute clustering
print "Compute structured hierarchical clustering..."
st = time.time()
ward = Ward(n_clusters=10).fit(X.T, adjacency_matrix)
label = ward.label_
print "Elaspsed time: ", time.time() - st
print "Size of label: ", label.shape[0]
print "Number of clusters: ", np.unique(label).shape[0]


###############################################################################
#Inverse the transformation to plot the results on an image
pl.figure()
axe = pl.subplot(1,2,1)
axe.imshow(np.reshape(label, mask.shape), interpolation="nearest",
                                                    cmap=pl.cm.jet)
axe = pl.subplot(1,2,2)
plot_dendrogram(axe, ward.parent_, ward.children_, ward.height_,
                active_nodes=ward.active_nodes_, )
pl.show()














