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
#categories = None

print "Loading 20 newsgroups dataset for categories:"
print categories

data_train = fetch_20newsgroups(subset='train', categories=categories,
                               shuffle=True, rng=42)

data_test = fetch_20newsgroups(subset='test', categories=categories,
                              shuffle=True, rng=42)

print "%d documents (training set)" % len(data_train.filenames)
print "%d documents (testing set)" % len(data_test.filenames)
print "%d categories" % len(data_train.target_names)
print

# split a training set and a test set
filenames_train, filenames_test = data_train.filenames, data_test.filenames
y_train, y_test = data_train.target, data_test.target

print "Extracting features from the training dataset using a sparse vectorizer"
t0 = time()
vectorizer = Vectorizer()
X_train = vectorizer.fit_transform((open(f).read() for f in filenames_train))
print "done in %fs" % (time() - t0)
print "n_samples: %d, n_features: %d" % X_train.shape
print

print "Extracting features from the test dataset using the same vectorizer"
t0 = time()
X_test = vectorizer.transform((open(f).read() for f in filenames_test))
print "done in %fs" % (time() - t0)
print "n_samples: %d, n_features: %d" % X_test.shape
print

# Build the affinity matrix using the cosnie similarity between documents
