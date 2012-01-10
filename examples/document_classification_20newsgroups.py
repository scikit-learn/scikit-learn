"""
======================================================
Classification of text documents using sparse features
======================================================

This is an example showing how the scikit-learn can be used to classify
documents by topics using a bag-of-words approach. This example uses
a scipy.sparse matrix to store the features instead of standard numpy arrays
and demos various classifiers that can efficiently handle sparse matrices.

The dataset used in this example is the 20 newsgroups dataset which will be
automatically downloaded and then cached.

You can adjust the number of categories by giving their names to the dataset
loader or setting them to None to get the 20 of them.

"""

# Author: Peter Prettenhofer <peter.prettenhofer@gmail.com>
#         Olivier Grisel <olivier.grisel@ensta.org>
#         Mathieu Blondel <mathieu@mblondel.org>
#         Lars Buitinck <L.J.Buitinck@uva.nl>
# License: Simplified BSD

import logging
import numpy as np
from operator import itemgetter
from optparse import OptionParser
import sys
from time import time

from sklearn.datasets import fetch_20newsgroups
from sklearn.feature_extraction.text import Vectorizer
from sklearn.feature_selection import SelectKBest, chi2
from sklearn.linear_model import RidgeClassifier
from sklearn.svm.sparse import LinearSVC
from sklearn.linear_model.sparse import SGDClassifier
from sklearn.naive_bayes import BernoulliNB, MultinomialNB
from sklearn.neighbors import KNeighborsClassifier
from sklearn.utils.extmath import density
from sklearn import metrics


# Display progress logs on stdout
logging.basicConfig(level=logging.INFO,
                    format='%(asctime)s %(levelname)s %(message)s')


# parse commandline arguments
op = OptionParser()
op.add_option("--report",
              action="store_true", dest="print_report",
              help="Print a detailed classification report.")
op.add_option("--chi2_select",
              action="store", type="int", dest="select_chi2",
              help="Select some number of features using a chi-squared test")
op.add_option("--confusion_matrix",
              action="store_true", dest="print_cm",
              help="Print the confusion matrix.")
op.add_option("--top10",
              action="store_true", dest="print_top10",
              help="Print ten most discriminative terms per class"
                   " for every classifier.")

(opts, args) = op.parse_args()
if len(args) > 0:
    op.error("this script takes no arguments.")
    sys.exit(1)

print __doc__
op.print_help()
print


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
print categories if categories else "all"

data_train = fetch_20newsgroups(subset='train', categories=categories,
                               shuffle=True, random_state=42)

data_test = fetch_20newsgroups(subset='test', categories=categories,
                              shuffle=True, random_state=42)
print 'data loaded'

categories = data_train.target_names    # for case categories == None

print "%d documents (training set)" % len(data_train.data)
print "%d documents (testing set)" % len(data_test.data)
print "%d categories" % len(categories)
print

# split a training set and a test set
y_train, y_test = data_train.target, data_test.target

print "Extracting features from the training dataset using a sparse vectorizer"
t0 = time()
vectorizer = Vectorizer()
X_train = vectorizer.fit_transform(data_train.data)
print "done in %fs" % (time() - t0)
print "n_samples: %d, n_features: %d" % X_train.shape
print

print "Extracting features from the test dataset using the same vectorizer"
t0 = time()
X_test = vectorizer.transform(data_test.data)
print "done in %fs" % (time() - t0)
print "n_samples: %d, n_features: %d" % X_test.shape
print

if opts.select_chi2:
    print ("Extracting %d best features by a chi-squared test" %
           opts.select_chi2)
    t0 = time()
    ch2 = SelectKBest(chi2, k=opts.select_chi2)
    X_train = ch2.fit_transform(X_train, y_train)
    X_test = ch2.transform(X_test)
    print "done in %fs" % (time() - t0)
    print

vocabulary = np.array([t for t, i in sorted(vectorizer.vocabulary.iteritems(),
                                            key=itemgetter(1))])


def trim(s):
    """Trim string to fit on terminal (assuming 80-column display)"""
    return s if len(s) <= 80 else s[:77] + "..."


###############################################################################
# Benchmark classifiers
def benchmark(clf):
    print 80 * '_'
    print "Training: "
    print clf
    t0 = time()
    clf.fit(X_train, y_train)
    train_time = time() - t0
    print "train time: %0.3fs" % train_time

    t0 = time()
    pred = clf.predict(X_test)
    test_time = time() - t0
    print "test time:  %0.3fs" % test_time

    score = metrics.f1_score(y_test, pred)
    print "f1-score:   %0.3f" % score

    if hasattr(clf, 'coef_'):
        print "dimensionality: %d" % clf.coef_.shape[1]
        print "density: %f" % density(clf.coef_)

        if opts.print_top10:
            print "top 10 keywords per class:"
            for i, category in enumerate(categories):
                top10 = np.argsort(clf.coef_[i])[-10:]
                print trim("%s: %s" % (category, " ".join(vocabulary[top10])))
        print

    if opts.print_report:
        print "classification report:"
        print metrics.classification_report(y_test, pred,
                                            target_names=categories)

    if opts.print_cm:
        print "confusion matrix:"
        print metrics.confusion_matrix(y_test, pred)

    print
    return score, train_time, test_time

for clf, name in ((RidgeClassifier(tol=1e-1), "Ridge Classifier"),
                  (KNeighborsClassifier(n_neighbors=10), "kNN")):
    print 80 * '='
    print name
    results = benchmark(clf)

for penalty in ["l2", "l1"]:
    print 80 * '='
    print "%s penalty" % penalty.upper()
    # Train Liblinear model
    liblinear_results = benchmark(LinearSVC(loss='l2', penalty=penalty, C=1000,
                                            dual=False, tol=1e-3))

    # Train SGD model
    sgd_results = benchmark(SGDClassifier(alpha=.0001, n_iter=50,
                                          penalty=penalty))

# Train SGD with Elastic Net penalty
print 80 * '='
print "Elastic-Net penalty"
sgd_results = benchmark(SGDClassifier(alpha=.0001, n_iter=50,
                                      penalty="elasticnet"))

# Train sparse Naive Bayes classifiers
print 80 * '='
print "Naive Bayes"
mnnb_results = benchmark(MultinomialNB(alpha=.01))
bnb_result = benchmark(BernoulliNB(alpha=.01))


class L1LinearSVC(LinearSVC):

    def fit(self, X, y):
        # The smaller C, the stronger the regularization.
        # The more regularization, the more sparsity.
        self.transformer_ = LinearSVC(C=1000, penalty="l1",
                                      dual=False, tol=1e-3)
        X = self.transformer_.fit_transform(X, y)
        return LinearSVC.fit(self, X, y)

    def predict(self, X):
        X = self.transformer_.transform(X)
        return LinearSVC.predict(self, X)

print 80 * '='
print "LinearSVC with L1-based feature selection"
l1linearsvc_results = benchmark(L1LinearSVC())
