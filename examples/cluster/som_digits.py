"""
=======================================================================
A demo of Self-Organising Map and KMeans on the handwritten digits data
=======================================================================

Comparing various SOM and Kmeans clustering on the handwritten digits data
with the pseudo_F index

"""
from __future__ import division
print __doc__

from time import time
import numpy as np

from sklearn.cluster import KMeans
from sklearn.cluster import SelfOrganizingMap
from sklearn.cluster import pseudo_F
from sklearn.datasets import load_digits
from sklearn.preprocessing import scale

np.random.seed(42)

# Load dataset

digits = load_digits()
data = scale(digits.data)
n_samples, n_features = data.shape
n_digits = len(np.unique(digits.target))

print "Digits dataset"
print "n_digits   : %d" % n_digits
print "n_features : %d" % n_features
print "n_samples  : %d" % n_samples
print

# Digits dataset clustering using Self-Organizing Map

print "Self-Organizing Map "
t0 = time()
grid_width = 4
som = SelfOrganizingMap(affinity=(grid_width, grid_width),
                        n_iterations=n_samples*5, learning_rate=1)
som.fit(data)
print "done in %0.3fs" % (time() - t0)
print

F = pseudo_F(data, som.labels_, som.centers_)
print 'pseudo_F %0.2f | %0.2f%%' % (F, 100 * (F / (1 + F)))
print

# Digits dataset clustering using Kmeans

print "KMeans "
t0 = time()
km = KMeans(init='k-means++', n_clusters=grid_width**2, n_init=10)
km.fit(data)
print "done in %0.3fs" % (time() - t0)
print

F = pseudo_F(data, km.labels_, km.cluster_centers_)
print 'pseudo_F %0.2f | %0.2f%%' % (F, 100 * (F / (1 + F)))
