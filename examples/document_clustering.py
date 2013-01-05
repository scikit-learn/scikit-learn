"""
=======================================
Clustering text documents using k-means
=======================================

This is an example showing how the scikit-learn can be used to cluster
documents by topics using a bag-of-words approach. This example uses
a scipy.sparse matrix to store the features instead of standard numpy arrays.

Two algorithms are demoed: ordinary k-means and its faster cousin minibatch
k-means.

"""

# Author: Peter Prettenhofer <peter.prettenhofer@gmail.com>
#         Lars Buitinck <L.J.Buitinck@uva.nl>
# License: Simplified BSD

from sklearn.datasets import fetch_20newsgroups
from sklearn.decomposition import LSA
from sklearn.feature_extraction.text import TfidfVectorizer
from sklearn import metrics

from sklearn.cluster import KMeans, MiniBatchKMeans

import logging
from optparse import OptionParser
import sys
from time import time

import numpy as np


# Display progress logs on stdout
logging.basicConfig(level=logging.INFO,
                    format='%(asctime)s %(levelname)s %(message)s')

# parse commandline arguments
op = OptionParser()
op.add_option("--lsa",
              dest="n_components", type="int",
              help="Preprocess documents with latent semantic analysis.")
op.add_option("--no-minibatch",
              action="store_false", dest="minibatch", default=True,
              help="Use ordinary k-means algorithm.")
op.add_option("--verbose",
              action="store_true", dest="verbose", default=False,
              help="Print progress reports inside k-means algorithm.")

print __doc__
op.print_help()

(opts, args) = op.parse_args()
if len(args) > 0:
    op.error("this script takes no arguments.")
    sys.exit(1)


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

labels = dataset.target
true_k = np.unique(labels).shape[0]

print "Extracting features from the training dataset using a sparse vectorizer"
t0 = time()
vectorizer = TfidfVectorizer(max_df=0.5, max_features=10000,
                             stop_words='english')
X = vectorizer.fit_transform(dataset.data)

print "done in %fs" % (time() - t0)
print

if opts.n_components:
    print "Performing dimensionality reduction using LSA"
    t0 = time()
    lsa = LSA(opts.n_components)
    X = lsa.fit_transform(X)

    print "done in %fs" % (time() - t0)
    print

print "n_samples: %d, n_features: %d" % X.shape
print


###############################################################################
# Do the actual clustering

if opts.minibatch:
    km = MiniBatchKMeans(n_clusters=true_k, init='k-means++', n_init=1,
                         init_size=1000, batch_size=1000, verbose=opts.verbose)
else:
    km = KMeans(n_clusters=true_k, init='random', max_iter=100, n_init=1,
                verbose=opts.verbose)

print "Clustering sparse data with\n    %r" % km
t0 = time()
km.fit(X)
print "done in %0.3fs" % (time() - t0)
print

print "Homogeneity: %0.3f" % metrics.homogeneity_score(labels, km.labels_)
print "Completeness: %0.3f" % metrics.completeness_score(labels, km.labels_)
print "V-measure: %0.3f" % metrics.v_measure_score(labels, km.labels_)
print "Adjusted Rand-Index: %.3f" % \
    metrics.adjusted_rand_score(labels, km.labels_)
print "Silhouette Coefficient: %0.3f" % metrics.silhouette_score(
    X, labels, sample_size=1000)

print
