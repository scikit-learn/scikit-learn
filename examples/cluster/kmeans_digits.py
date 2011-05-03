"""
===========================================================
A demo of K-Means clustering on the handwritten digits data
===========================================================

Comparing various initialization strategies in terms of runtime and quality of
the results.

TODO: explode the ouput of the cluster labeling and digits.target groundtruth
as categorical boolean arrays of shape (n_sample, n_unique_labels) and measure
the Pearson correlation as an additional measure of the clustering quality.
"""
print __doc__

from time import time
import numpy as np

from scikits.learn.cluster import KMeans
from scikits.learn.datasets import load_digits
from scikits.learn.decomposition import PCA
from scikits.learn.preprocessing import scale

np.random.seed(42)

digits = load_digits()
data = scale(digits.data)

n_samples, n_features = data.shape
n_digits = len(np.unique(digits.target))

print "n_digits: %d" % n_digits
print "n_features: %d" % n_features
print "n_samples: %d" % n_samples
print

print "Raw k-means with k-means++ init..."
t0 = time()
km = KMeans(init='k-means++', k=n_digits, n_init=10).fit(data)
print "done in %0.3fs" % (time() - t0)
print "inertia: %f" % km.inertia_
print

print "Raw k-means with random centroid init..."
t0 = time()
km = KMeans(init='random', k=n_digits, n_init=10).fit(data)
print "done in %0.3fs" % (time() - t0)
print "inertia: %f" % km.inertia_
print

print "Raw k-means with PCA-based centroid init..."
# in this case the seeding of the centers is deterministic, hence we run the
# kmeans algorithm only once with n_init=1
t0 = time()
pca = PCA(n_components=n_digits).fit(data)
km = KMeans(init=pca.components_, k=n_digits, n_init=1).fit(data)
print "done in %0.3fs" % (time() - t0)
print "inertia: %f" % km.inertia_
print

