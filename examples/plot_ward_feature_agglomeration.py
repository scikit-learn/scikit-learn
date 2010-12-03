"""
===========================================================
A demo of feature agglomeration - structured ward
===========================================================

Author : Vincent Michel, 2010

"""
print __doc__


import numpy as np
import pylab as pl
from scikits.learn.feature_agglomeration import WardAgglomeration
from scikits.learn.cluster import plot_dendrogram
import scipy.linalg as sl
import time
from scikits.learn.feature_extraction.image import img_to_graph
from scikits.learn.linear_model import BayesianRidge


###############################################################################
# Generate data
n_samples = 200
size = 100
roi_size = 10
snr = 5.
np.random.seed(0)
mask = np.ones([size, size], dtype=np.bool)
w = np.zeros((size, size))
w[0:roi_size, 0:roi_size] = -1.
w[-roi_size:, -roi_size:] = 1.
X = np.random.randn(n_samples, size * size)
w = w.ravel()
y_no_noise = np.array(np.dot(X, w))
orig_noise = np.random.randn(y_no_noise.shape[0])
noise_coef = (sl.norm(y_no_noise, 2) / np.exp(snr / 20.))\
             / sl.norm(orig_noise, 2)
y = y_no_noise + noise_coef * orig_noise
X -= X.mean(axis=-1)[:, np.newaxis]
X /= X.std(axis=-1)[:, np.newaxis]

n_clusters = 30

###############################################################################
# Feature agglomeration
print "Compute structured hierarchical clustering..."
st = time.time()
adjacency_matrix = img_to_graph(mask,mask)
feat_agglo = WardAgglomeration(n_clusters)
feat_agglo.fit(X, adjacency_matrix)
label = feat_agglo.label_
Xred = feat_agglo.transform(X, pooling_func=np.mean)
print "Elaspsed time: ", time.time() - st
print "Size of label: ", label.shape[0]
print "Number of clusters: ", np.unique(label).shape[0]
print "Initial data shape: ", X.shape
print "Reduced data shape: ", Xred.shape

###############################################################################
# Compute the weights of a Bayesian Ridge
raw_coef_ = BayesianRidge().fit(Xred, y).coef_
coef_ = np.reshape(feat_agglo.inverse_transform(raw_coef_), mask.shape[:2])

###############################################################################
#Inverse the transformation to plot the results on an image
pl.figure()
axe = pl.subplot(1,3,1)
axe.imshow(np.reshape(w, mask.shape[:2]),
interpolation="nearest", cmap=pl.cm.RdBu_r)
axe.set_title("True weights")
axe = pl.subplot(1,3,2)
axe.imshow(coef_, interpolation="nearest", cmap=pl.cm.RdBu_r)
axe.set_title("Estimated weights")
axe = pl.subplot(1,3,3)
ward = feat_agglo.clustering
plot_dendrogram(axe, ward.parent_, ward.children_, ward.height_,
            active_nodes=ward.active_nodes_,weights_nodes=raw_coef_,
            cmap_nodes=pl.cm.RdBu_r)
pl.show()












