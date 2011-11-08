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

from sklearn.datasets import fetch_20newsgroups
from sklearn.feature_extraction.text import Vectorizer
from sklearn.cluster import MiniBatchKMeans
from sklearn.cluster import power_iteration_clustering
from sklearn.cluster import spectral_clustering
from sklearn.metrics.pairwise import cosine_similarity
from sklearn.metrics import homogeneity_completeness_v_measure
from sklearn.metrics import adjusted_rand_score

# Display progress logs on stdout
logging.basicConfig(level=logging.INFO,
                    format='%(asctime)s %(levelname)s %(message)s')

random_state = None

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
                                shuffle=True, random_state=0)

print "%d documents (training set)" % len(data_train.filenames)
print "%d categories" % len(data_train.target_names)
print

from sklearn.externals import joblib
memory = joblib.Memory('.')

@memory.cache
def extract(documents):
    print "Extracting features from the training dataset using a vectorizer"
    t0 = time()
    vectorizer = Vectorizer()
    X = vectorizer.fit_transform(documents)
    print "done in %fs" % (time() - t0)
    print
    return X

X_train = extract(data_train.data)
labels_true = data_train.target

# Smaller dataset to reproduce conditions similar to the PIC paper
X_train = X_train[:800]
labels_true = labels_true[:800]
print "n_samples: %d, n_features: %d" % X_train.shape

# Build the affinity matrix using the cosine similarity between documents

print "Extracting the pairwise cosine similarity between documents"
t0 = time()
affinity = cosine_similarity(X_train, X_train)
print "done in %fs" % (time() - t0)
print

n_clusters = len(data_train.target_names)

def report(labels_true, labels_pred, duration):
    """Print clustering report on stdout"""
    h, c, v = homogeneity_completeness_v_measure(labels_true, labels_pred)
    print "Homogeneity: %0.3f" % h
    print "Completeness: %0.3f" % c
    print "V-Measure: %0.3f" % v
    print "ARI: %0.3f" % adjusted_rand_score(labels_true, labels_pred)
    print "Duration: %0.3fs" % duration
    print


print "Clustering the documents using the MiniBatchKMeans"
t0 = time()
labels_pred_pic = MiniBatchKMeans(
    k=n_clusters, verbose=2, random_state=random_state).fit(X_train).labels_
duration = time() - t0
report(labels_true, labels_pred_pic, duration)

print "Clustering the documents using the Power Iteration Clustering method"
t0 = time()
labels_pred_pic = power_iteration_clustering(
    affinity, k=n_clusters, verbose=2, n_vectors=1, tol=1e-6,
    random_state=random_state)
duration = time() - t0
report(labels_true, labels_pred_pic, duration)

print "Clustering the documents using the Spectral Clustering method"
t0 = time()
labels_pred_spectral = spectral_clustering(
    affinity, k=n_clusters, mode='arpack', random_state=random_state)
duration = time() - t0
report(labels_true, labels_pred_spectral, duration)
#
#print "Agreement between PIC and Spectral labelings"
#report(labels_pred_pic, labels_pred_spectral, 0.0)
