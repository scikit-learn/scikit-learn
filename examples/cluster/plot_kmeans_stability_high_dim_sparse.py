"""
===================================================================
Empirical evaluation of the impact of k-means initialization (text)
===================================================================

Evaluate the ability of k-means initializations strategies to make
the algorithm convergence robust as measured by the relative standard
deviation of the inertia of the clustering (i.e. the sum of distances
to the nearest cluster center).

The dataset used for evaluation 4 categories of the 20 newsgroups
dataset.

"""
print __doc__

# Author: Olivier Grisel <olivier.grisel@ensta.org>
# License: Simplified BSD

import numpy as np
import pylab as pl
import matplotlib.cm as cm
from time import time

from sklearn.cluster import MiniBatchKMeans
from sklearn.datasets import fetch_20newsgroups
from sklearn.feature_extraction.text import Vectorizer

random_state = np.random.RandomState(0)

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
vectorizer = Vectorizer(max_df=0.95, max_features=10000)
X = vectorizer.fit_transform(dataset.data)

print "done in %fs" % (time() - t0)
print "n_samples: %d, n_features: %d" % X.shape
print

# Number of run (with randomly generated dataset) for each strategy so as
# to be able to compute an estimate of the standard deviation
n_runs = 5

# k-means models can do several random inits so as to be able to trade
# CPU time for convergence robustness
n_init_range = np.array([1, 5, 10])

fig = pl.figure()
plots = []
legends = []

cases = [
    ('k-means++', {'max_no_improvement': 5, 'n_reinit': 0}),
    ('k-means++', {'max_no_improvement': 5, 'n_reinit': 2}),
    ('random', {'max_no_improvement': 5, 'n_reinit': 0}),
    ('random', {'max_no_improvement': 5, 'n_reinit': 2}),
]

for init, params in cases:
    print "Evaluation of batch MiniBatchKMeans with %s init" % init
    inertia = np.empty((len(n_init_range), n_runs))

    for run_id in range(n_runs):
        for i, n_init in enumerate(n_init_range):
            km = MiniBatchKMeans(k=true_k,
                                 init=init,
                                 random_state=run_id,
                                 n_init=n_init,
                                 **params).fit(X)
            inertia[i, run_id] = km.inertia_
            print "Inertia for n_init=%02d, run_id=%d: %0.3f" % (
                n_init, run_id, km.inertia_)

    plots.append(
        pl.errorbar(n_init_range, inertia.mean(axis=1), inertia.std(axis=1)))
    n_reinit = params.get('n_reinit')
    if n_reinit is not None:
        legends.append("MiniBatchKMeans with %s init and %d reinit" % (
            init, n_reinit))

pl.xlabel('n_init')
pl.ylabel('inertia')
pl.legend(plots, legends)
pl.title("Mean inertia for various k-means init accross %d runs" % n_runs)
pl.show()
