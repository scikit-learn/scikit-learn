"""
==================================================
Clustering of text documents using sparse features
==================================================

This is an example showing how the scikit-learn can be used to cluster
documents by topics using a bag-of-words approach. This example uses
a scipy.sparse matrix to store the features instead of standard numpy arrays.

The dataset used in this example is the 20 newsgroups dataset which will be
automatically downloaded and then cached.

You can adjust the number of categories by giving there name to the dataset
loader or setting them to None to get the 20 of them.

This example demos various linear classifiers with different training
strategies.

"""
print __doc__

# Author: Olivier Grisel <olivier.grisel@ensta.org>
# License: Simplified BSD

from time import time
import logging
import os
import sys

from scikits.learn.datasets import fetch_20newsgroups
from scikits.learn.feature_extraction.text import Vectorizer
from scikits.learn.cluster import power_iteration_clustering
from scikits.learn.cluster import spectral_clustering
from scikits.learn.metrics.pairwise import cosine_similarity

# Display progress logs on stdout
logging.basicConfig(level=logging.INFO,
                    format='%(asctime)s %(levelname)s %(message)s')

################################################################################
# Load some categories from the training set
categories = [
    'misc.forsale',
    'soc.religion.christian',
    'talk.politics.guns',
    'rec.sport.baseball',
]
# Uncomment the following to do the analysis on all the categories
#categories = None

print "Loading 20 newsgroups dataset for categories:"
print categories

data_train = fetch_20newsgroups(subset='train', categories=categories,
                                shuffle=True, rng=42)

print "%d documents (training set)" % len(data_train.filenames)
print "%d categories" % len(data_train.target_names)
print

filenames_train = data_train.filenames
y_train = data_train.target

print "Extracting features from the training dataset using a sparse vectorizer"
t0 = time()
vectorizer = Vectorizer()
X_train = vectorizer.fit_transform((open(f).read() for f in filenames_train))
X_train = X_train[:800] # try to reproduce the settings of the PIC paper
print "done in %fs" % (time() - t0)
print "n_samples: %d, n_features: %d" % X_train.shape
print

# Build the affinity matrix using the cosine similarity between documents

print "Extracting the pairwise cosine similarity between documents"
t0 = time()
affinity = cosine_similarity(X_train, X_train)
print "done in %fs" % (time() - t0)
print

n_clusters = len(data_train.target_names)

print "Clustering the documents using the Power Iteration Clustering method"
t0 = time()
labels = power_iteration_clustering(affinity, k=n_clusters, n_vectors=5)
print "done in %fs" % (time() - t0)
for i in range(20):
    print "Doc %d from category '%s' assigned to cluster %d" % (
        i, data_train.target_names[data_train.target[i]], labels[i])
print


print "Clustering the documents using the Spectral Clustering method"
t0 = time()
labels = spectral_clustering(affinity, k=n_clusters, mode='arpack')
print "done in %fs" % (time() - t0)
for i in range(20):
    print "Doc %d from category '%s' assigned to cluster %d" % (
        i, data_train.target_names[data_train.target[i]], labels[i])
print
