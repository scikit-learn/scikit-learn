"""
=======================================================
Comparison of PCA and CCIPCA on a 2D dataset
=======================================================

The dataset is a cluster of 2D points that is roughly
ellipitical.

Principal Component Analysis (PCA) applied to this data identifies the
primary and secondary axis of the ellipse.

Candid Covariance-free Principal Component Analysis CCIPCA identifies
the same axes but processed the data in a stream rather than in a batch.
"""

import pylab as pl
import numpy as np
import math

from sklearn.decomposition import PCA
from sklearn.decomposition import CCIPCA

# 2D point set in the shape of an ellipse
X = np.array([10*np.random.randn(1000), np.random.randn(1000)]).transpose()
rot = np.array( [[math.cos(math.radians(135)), -math.sin(math.radians(135))],
                [math.sin(math.radians(135)),  math.cos(math.radians(135))]] )
X = np.dot(X,rot)

# PCA
pca = PCA(n_components=2)
pca.fit(X)
print 'pca'
print pca.components_
print pca.explained_variance_ratio_

# CCIPCA
ccipca = CCIPCA(n_components=2)   
ccipca.fit(X)
print 'ccipca'
print ccipca.components_
print ccipca.explained_variance_ratio_

# Plot PCA
pl.figure()
# the point data set
pl.scatter( X[:,0], X[:,1] )
# the first component
pl.arrow(pca.mean_[0], pca.mean_[1], 
            25*pca.components_[0,0], 25*pca.components_[0,1], 
           shape='full', lw=3,length_includes_head=True, head_width=1.25, color = 'r')
# the second component
pl.arrow(pca.mean_[0], pca.mean_[1], 
        -5*pca.components_[1,0], -5*pca.components_[1,1], 
        shape='full', lw=3,length_includes_head=True, head_width=1.25, color = 'r')
pl.title('PCA of 2D ellipse')

# Plot CCIPCA
pl.figure()
# the point data set
pl.scatter( X[:,0], X[:,1] )
# the first component
pl.arrow(ccipca.mean_[0], ccipca.mean_[1], 
            25*ccipca.components_[0,0], 25*ccipca.components_[0,1], 
           shape='full', lw=3,length_includes_head=True, head_width=1.25, color = 'r')
# the second component
pl.arrow(ccipca.mean_[0], ccipca.mean_[1], 
        -5*ccipca.components_[1,0], -5*ccipca.components_[1,1], 
        shape='full', lw=3,length_includes_head=True, head_width=1.25, color = 'r')
pl.title('CCIPCA of 2D ellipse')

pl.show()
