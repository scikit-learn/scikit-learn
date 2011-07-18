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
from scikits.learn.cluster import randindex

from scikits.learn.preprocessing import Normalizer


# Display progress logs on stdout
logging.basicConfig(level=logging.INFO,
                    format='%(asctime)s %(levelname)s %(message)s')

###############################################################################
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
vectorizer = Vectorizer(max_features=10000)
X = vectorizer.fit_transform((open(f).read() for f in filenames))

X = Normalizer(norm="l2", copy=False).transform(X)

print "done in %fs" % (time() - t0)
print "n_samples: %d, n_features: %d" % X.shape
print


###############################################################################
# Now sparse MiniBatchKmeans

print "_" * 80

mbkm = MiniBatchKMeans(init="random", k=true_k, max_iter=10, random_state=13,
                       chunk_size=1000, tol=0.0, n_init=1)

print "Clustering sparse data with %s" % str(mbkm)
print

t0 = time()
mbkm.fit(X)
print "done in %0.3fs" % (time() - t0)

ri = randindex(labels, mbkm.labels_)
vmeasure = metrics.v_measure_score(labels, mbkm.labels_)
print "Homogeneity: %0.3f" % metrics.homogeneity_score(labels, mbkm.labels_)
print "Completeness: %0.3f" % metrics.completeness_score(labels, mbkm.labels_)
print "V-measure: %0.3f" % vmeasure
print "Rand-Index: %.3f" % ri
print
