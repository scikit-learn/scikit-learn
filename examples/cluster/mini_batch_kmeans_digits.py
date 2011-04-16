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
from scikits.learn.preprocessing import scale

np.random.seed(42)

digits = load_digits()
data = scale(digits.data)

n_samples, n_features = data.shape
n_digits = len(np.unique(digits.target))

def mini_batch(X, batch_size, k = 10, iterations = 300):
    km = MiniBatchKMeans(init='k-means++',
                         k=k)

    for i in xrange(iterations):
        j = i * batch_size % len(data)
        sample = X[j:j+batch_size]
        km.partial_fit(sample)

    # That last calculation is done only to get the inertia, to be able to
    # compare to the previous calculation
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
km = mini_batch(data, 600, k=n_digits)
print "done in %0.3fs" % (time() - t0)
print "inertia: %f" % km.inertia_
print


print "Mini batch k-means with k-means++ init, chunk 300..."
t0 = time()
km = mini_batch(data, 300, k=n_digits)
print "done in %0.3fs" % (time() - t0)
print "inertia: %f" % km.inertia_
print

print "Mini batch k-means with k-means++ init, chunk 100..."
t0 = time()
km = mini_batch(data, 100, k=n_digits)
print "done in %0.3fs" % (time() - t0)
print "inertia: %f" % km.inertia_
print

print "Mini batch k-means with k-means++ init, chunk 20..."
t0 = time()
km = mini_batch(data, 20, k=n_digits)
print "done in %0.3fs" % (time() - t0)
print "inertia: %f" % km.inertia_
print


