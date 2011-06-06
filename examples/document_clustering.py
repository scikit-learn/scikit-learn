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
from scikits.learn.cluster import sparse
from scikits.learn.cluster.sparse import randindex

from scikits.learn.preprocessing.sparse import LengthNormalizer


## def randindex(labels_true, labels_pred):
##     count = np.zeros((2, 2), dtype=np.int32)
##     assert len(labels_true) == len(labels_pred)

##     n = len(labels_pred)
##     if n < 2:
##         raise ValueError("number of samples must be at least 2.")

##     for i in xrange(n):
##         for j in xrange(i, n):
##             if labels_true[i] == labels_true[j]:
##                 if labels_pred[i] == labels_pred[j]:
##                     count[1, 1] += 1
##                 else:
##                     count[1, 0] += 1
##             else:
##                 if labels_pred[i] == labels_pred[j]:
##                     count[0, 1] += 1
##                 else:
##                     count[0, 0] += 1

##     return float(count[0, 0] + count[1, 1]) / float((n * (n - 1)) / 2.0)



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

print "%d documents (training set)" % len(data_train.filenames)
print "%d categories" % len(data_train.target_names)
print

# split a training set and a test set
filenames_train = data_train.filenames
labels = data_train.target
true_k = np.unique(labels).shape[0]

print "Extracting features from the training dataset using a sparse vectorizer"
t0 = time()
vectorizer = Vectorizer(max_features=10000)
X = vectorizer.fit_transform((open(f).read() for f in filenames_train))

X = LengthNormalizer().transform(X)

print "done in %fs" % (time() - t0)
print "n_samples: %d, n_features: %d" % X.shape
print

chunk_size = 300

print "_" * 80

mbkm = sparse.MiniBatchKMeans(k=true_k, n_iter=100, random_state=13,
                              chunk_size=chunk_size)

print "Clustering data with %s" % str(mbkm)
print

t0 = time()

mbkm.fit(X)

print "done in %0.3fs" % (time() - t0)
print "Homogeneity: %0.3f" % metrics.homogeneity_score(labels, mbkm.labels_)
print "Completeness: %0.3f" % metrics.completeness_score(labels, mbkm.labels_)
print "V-measure: %0.3f" % metrics.v_measure_score(labels, mbkm.labels_)
print "Rand-Index: %.3f" % randindex(labels, mbkm.labels_)
print

## ################################################################################
## # Now dense Kmeans

## print "_" * 80

## km = KMeans(k=true_k, init="k-means++", n_init=1, random_state=13, tol=0.001,
##             copy_x=False)

## X = X.toarray()

## print "Clustering data with %s" % str(km)
## print

## t0 = time()

## km.fit(X)

## print "done in %0.3fs" % (time() - t0)
## print "Homogeneity: %0.3f" % metrics.homogeneity_score(labels, km.labels_)
## print "Completeness: %0.3f" % metrics.completeness_score(labels, km.labels_)
## print "V-measure: %0.3f" % metrics.v_measure_score(labels, km.labels_)
## print "Rand-Index: %.3f" % randindex(labels, km.labels_)
## print

## ################################################################################
## # Now dense MiniBatchKmeans

## print "_" * 80

## mbkm = MiniBatchKMeans(k=true_k, init="random", n_init=1, max_iter=100,
##                        random_state=13, chunk_size=chunk_size)

## print "Clustering data with %s" % str(mbkm)
## print

## t0 = time()

## mbkm.fit(X)

## print "done in %0.3fs" % (time() - t0)
## print "Homogeneity: %0.3f" % metrics.homogeneity_score(labels, mbkm.labels_)
## print "Completeness: %0.3f" % metrics.completeness_score(labels, mbkm.labels_)
## print "V-measure: %0.3f" % metrics.v_measure_score(labels, mbkm.labels_)
## print "Rand-Index: %.3f" % randindex(labels, mbkm.labels_)
## print
