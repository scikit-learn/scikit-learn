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

from sklearn.datasets import fetch_20newsgroups
from sklearn.feature_extraction.text import Vectorizer
from sklearn import metrics

from sklearn.cluster import MiniBatchKMeans



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
#categories = None

print "Loading 20 newsgroups dataset for categories:"
print categories

dataset = fetch_20newsgroups(subset='all', categories=categories,
                             shuffle=True, random_state=42)

print "%d documents" % len(dataset.data)
print "%d categories" % len(dataset.target_names)
print

# split a training set and a test set
labels = dataset.target
true_k = np.unique(labels).shape[0]

print "Extracting features from the training dataset using a sparse vectorizer"
t0 = time()
vectorizer = Vectorizer(max_df=0.95, max_features=10000)
X = vectorizer.fit_transform(dataset.data)

print "done in %fs" % (time() - t0)
print "n_samples: %d, n_features: %d" % X.shape
print

###############################################################################
# Now sparse MiniBatchKmeans

mbkm = MiniBatchKMeans(init="random", k=true_k, max_iter=10, chunk_size=1000,
                       verbose=0)
print "Clustering sparse data with %s" % mbkm
t0 = time()
mbkm.fit(X)
print "done in %0.3fs" % (time() - t0)
print

print "Homogeneity: %0.3f" % metrics.homogeneity_score(labels, mbkm.labels_)
print "Completeness: %0.3f" % metrics.completeness_score(labels, mbkm.labels_)
print "V-measure: %0.3f" % metrics.v_measure_score(labels, mbkm.labels_)
print "Adjusted Rand-Index: %.3f" % \
    metrics.adjusted_rand_score(labels, mbkm.labels_)
print "Silhouette Coefficient: %0.3f" % metrics.silhouette_score(
    X, labels, sample_size=1000)

print
