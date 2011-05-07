"""
========================================================================
A comparison of minibatch K-Means and K-Means on handwritten digits data
========================================================================

Comparing the minibatch K-Means with the normal K-Means algorithm in terms of
runtime and quality of the results.
"""
print __doc__

from time import time
import numpy as np

from scikits.learn.cluster import MiniBatchKMeans, KMeans
from scikits.learn.datasets import load_digits
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

print "Mini batch k-means with k-means++ init, chunk 600..."
t0 = time()
km = MiniBatchKMeans(init='k-means++',
                     k=n_digits,
                     chunk_size=600).fit(data)
print "done in %0.3fs" % (time() - t0)
print "inertia: %f" % km.inertia_
print


print "Mini batch k-means with k-means++ init, chunk 300..."
t0 = time()
km = MiniBatchKMeans(init='k-means++',
                     k=n_digits,
                     chunk_size=300).fit(data)
print "done in %0.3fs" % (time() - t0)
print "inertia: %f" % km.inertia_
print

print "Mini batch k-means with k-means++ init, chunk 100..."
t0 = time()
km = MiniBatchKMeans(init='k-means++',
                     k=n_digits,
                     chunk_size=100).fit(data)
print "done in %0.3fs" % (time() - t0)
print "inertia: %f" % km.inertia_
print

print "Mini batch k-means with k-means++ init, chunk 20..."
t0 = time()
km = MiniBatchKMeans(init='k-means++',
                     k=n_digits,
                     chunk_size=20).fit(data)
print "done in %0.3fs" % (time() - t0)
print "inertia: %f" % km.inertia_
print
