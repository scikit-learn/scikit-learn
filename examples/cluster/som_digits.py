"""
=======================================================================
A demo of Self-Organising Map and KMeans on the handwritten digits data
=======================================================================

XXX : Should add text to describe what the example and what to expect
from the output. Would it be possible to plot something?
"""
from __future__ import division
print __doc__

from time import time
import numpy as np

from scikits.learn.cluster import KMeans
from scikits.learn.cluster import SelfOrganizingMap
from scikits.learn.cluster import calinski_index

from scikits.learn.datasets import load_digits
from scikits.learn.preprocessing import scale


def display(labels, digits, n_clusters):
    # XXX : n_clusters unused
    r = {0: [], 1: [], 2: [], 3: [], 4: [], 5: [], 6: [], 7: [], 8: [], 9: []}
    for i, v in enumerate(labels):
        r[digits.target[i]].append(v)

    for k, v in r.items():
        s = set(v)
        print 'target %i | nb cluster %i |' % (k, len(s)), s

np.random.seed(42)

digits = load_digits()
data = scale(digits.data)

n_samples, n_features = data.shape
n_digits = len(np.unique(digits.target))

print "n_digits: %d" % n_digits
print "n_features: %d" % n_features
print "n_samples: %d" % n_samples
print

print "Self-Organizing Map "
t0 = time()
grid_width = 4
som = SelfOrganizingMap(size=grid_width, n_iterations=n_samples*5,
                        learning_rate=1)
som.fit(data)
print "done in %0.3fs" % (time() - t0)
print

display(som.labels_, digits, grid_width**2)
C = calinski_index(data, som.labels_, som.neurons_)
print 'calinski index %0.2f | %0.2f%%' % (C, 100 * (C / (1 + C)))
print

print "KMeans "
t0 = time()
km = KMeans(init='k-means++', k=grid_width**2, n_init=10)
km.fit(data)
print "done in %0.3fs" % (time() - t0)
print

display(km.labels_, digits, n_digits)
C = calinski_index(data, km.labels_, km.cluster_centers_)
print 'calinski index %0.2f | %0.2f%%' % (C, 100 * (C / (1 + C)))
