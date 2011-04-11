"""
===========================================================================
A comparison of batch K-Means and normal K-Means on handwritten digits data
===========================================================================

Comparing the batch K-Means with the normal K-Means algorithm in terms of
runtime and quality of the results.
"""
print __doc__

from time import time
import numpy as np

from scikits.learn.cluster import MiniBatchKMeans, KMeans
from scikits.learn.datasets import load_digits
from scikits.learn.pca import PCA
from scikits.learn.preprocessing import scale

np.random.seed(42)

digits = load_digits()
data = scale(digits.data)

n_samples, n_features = data.shape
n_digits = len(np.unique(digits.target))

def mini_batch(X, chunk, iterations):
    km = MiniBatchKMeans(init='k-means++',
                         k=n_digits,
                         n_init=10)

    for i in xrange(iterations):
        j = i * chunk % len(data)
        sample = data[j:j+chunk]
        km.partial_fit(sample)

    # That last calculation is done only to get the inertia, to be able to compare
    # to the previous calculation
    km.partial_fit(data)
    return km


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
km = mini_batch(data, 600, 300)
print "done in %0.3fs" % (time() - t0)
print "inertia: %f" % km.inertia_
print


print "Mini batch k-means with k-means++ init, chunk 300..."
t0 = time()
km = mini_batch(data, 600, 300)
print "done in %0.3fs" % (time() - t0)
print "inertia: %f" % km.inertia_
print

print "Mini batch k-means with k-means++ init, chunk 150..."
t0 = time()
km = mini_batch(data, 600, 300)
print "done in %0.3fs" % (time() - t0)
print "inertia: %f" % km.inertia_
print

print "Mini batch k-means with k-means++ init, chunk 100..."
t0 = time()
km = mini_batch(data, 600, 300)
print "done in %0.3fs" % (time() - t0)
print "inertia: %f" % km.inertia_
print


