"""
===============================================
Clustering text documents using MiniBatchKmeans
===============================================

This is an example showing how the scikit-learn can be used to cluster
documents by topics using a bag-of-words approach. This example uses
a scipy.sparse matrix to store the features instead of standard numpy arrays.

"""
print __doc__

# Author: Peter Prettenhofer <peter.prettenhofer@gmail.com>
# License: Simplified BSD

from time import time
import logging
import numpy as np

from scikits.learn.datasets import fetch_20newsgroups
from scikits.learn.feature_extraction.text import Vectorizer
from scikits.learn import metrics

from scikits.learn.cluster import MiniBatchKMeans
from scikits.learn.cluster import KMeans
from scikits.learn.cluster import SpectralClustering
from scikits.learn.cluster import sparse
from scikits.learn.cluster.sparse import randindex

from scikits.learn.preprocessing import Normalizer


# Display progress logs on stdout
logging.basicConfig(level=logging.INFO,
                    format='%(asctime)s %(levelname)s %(message)s')

################################################################################
# Load some categories from the training set
categories = [
    'alt.atheism',
    'talk.religion.misc',
    'comp.graphics',
    'sci.space',
]
# Uncomment the following to do the analysis on all the categories
# categories = None

print "Loading 20 newsgroups dataset for categories:"
print categories

data_train = fetch_20newsgroups(subset='train', categories=categories,
                               shuffle=True, random_state=42)
data_test = fetch_20newsgroups(subset='test', categories=categories,
                               shuffle=True, random_state=42)

filenames = np.concatenate((data_train.filenames, data_test.filenames))
target_names = set(data_train.target_names + data_test.target_names)

print "%d documents" % len(filenames)
print "%d categories" % len(target_names)
print

# split a training set and a test set
labels = np.concatenate((data_train.target, data_test.target))
true_k = np.unique(labels).shape[0]

print "Extracting features from the training dataset using a sparse vectorizer"
t0 = time()
vectorizer = Vectorizer(max_features=1000)
X = vectorizer.fit_transform((open(f).read() for f in filenames))

X = Normalizer(norm="l2", copy=False).transform(X)

print "done in %fs" % (time() - t0)
print "n_samples: %d, n_features: %d" % X.shape
print

chunk_size = 500


################################################################################
# Now sparse MiniBatchKmeans

print "_" * 80

mbkm = sparse.MiniBatchKMeans(k=true_k, max_iter=100, random_state=13,
                              chunk_size=chunk_size, tol=0.0)

print "Clustering data with %s" % str(mbkm)
print

t0 = time()

mbkm.fit(X)
ri = randindex(labels, mbkm.labels_)
vmeasure = metrics.v_measure_score(labels, mbkm.labels_)
print "done in %0.3fs" % (time() - t0)
print "Homogeneity: %0.3f" % metrics.homogeneity_score(labels, mbkm.labels_)
print "Completeness: %0.3f" % metrics.completeness_score(labels, mbkm.labels_)
print "V-measure: %0.3f" % vmeasure
print "Rand-Index: %.3f" % ri
print "center norms", np.sum(mbkm.cluster_centers_ ** 2.0, axis=1)
print


################################################################################
# Now dense MiniBatchKmeans

print "_" * 80

mbkm = MiniBatchKMeans(init="random", k=true_k, max_iter=100, random_state=13,
                       chunk_size=chunk_size, tol=0.0, n_init=1)

print "Clustering data with %s" % str(mbkm)
print

t0 = time()

mbkm.fit(X.toarray())
ri = randindex(labels, mbkm.labels_)
vmeasure = metrics.v_measure_score(labels, mbkm.labels_)
print "done in %0.3fs" % (time() - t0)
print "Homogeneity: %0.3f" % metrics.homogeneity_score(labels, mbkm.labels_)
print "Completeness: %0.3f" % metrics.completeness_score(labels, mbkm.labels_)
print "V-measure: %0.3f" % vmeasure
print "Rand-Index: %.3f" % ri
print "center norms", np.sum(mbkm.cluster_centers_ ** 2.0, axis=1)
print



## ################################################################################
## # Now sectral clustering

## print "_" * 80
## sc = SpectralClustering(k=true_k, mode='amg', random_state=13)

## t0 = time()
## A = (X * X.T).toarray()
## print "computed A in %0.3fs" % (time() - t0)

## print "Clustering data with %s" % str(sc)
## print

## t0 = time()

## sc.fit(A)

## print "done in %0.3fs" % (time() - t0)
## print "Homogeneity: %0.3f" % metrics.homogeneity_score(labels, sc.labels_)
## print "Completeness: %0.3f" % metrics.completeness_score(labels, sc.labels_)
## print "V-measure: %0.3f" % metrics.v_measure_score(labels, sc.labels_)
## print "Rand-Index: %.3f" % randindex(labels, sc.labels_)
## print 
