"""
===============================================
Semisupervised classification of text documents
===============================================

This variation on the document classification theme (see
document_classification_20newsgroups.py) showcases semisupervised learning:
classification with training on partially unlabeled data.

The dataset used in this example is the 20 newsgroups dataset which will be
automatically downloaded and then cached; this set is labeled, but the
labels from a random part will be removed.

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

from sklearn.cross_validation import StratifiedKFold
from sklearn.datasets import fetch_20newsgroups
from sklearn.feature_extraction.text import Vectorizer
from sklearn.naive_bayes import BernoulliNB, SemisupervisedNB, MultinomialNB
from sklearn import metrics


# Display progress logs on stdout
logging.basicConfig(level=logging.INFO,
                    format='%(asctime)s %(levelname)s %(message)s')


# parse commandline arguments
op = OptionParser()
op.add_option("--confusion_matrix",
              action="store_true", dest="print_cm",
              help="Print the confusion matrix.")
op.add_option("--labeled",
              action="store", type="float", dest="labeled_fraction",
              help="Fraction of labels to retain (roughly).")
op.add_option("--report",
              action="store_true", dest="print_report",
              help="Print a detailed classification report.")
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


def split_indices(y, fraction):
    """Random stratified split of indices into y
    
    Returns (unlabeled, labeled)
    """
    k = int(round(1 / fraction))
    folds = list(StratifiedKFold(y, k))
    return folds[rng.randint(k)]


def trim(s):
    """Trim string to fit on terminal (assuming 80-column display)"""
    return s if len(s) <= 80 else s[:77] + "..."


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

rng = np.random.RandomState(42)

data_train = fetch_20newsgroups(subset='train', categories=categories,
                               shuffle=True, random_state=rng)

data_test = fetch_20newsgroups(subset='test', categories=categories,
                              shuffle=True, random_state=rng)
print 'data loaded'

categories = data_train.target_names    # for case categories == None

print "%d documents (training set)" % len(data_train.data)
print "%d documents (testing set)" % len(data_test.data)
print "%d categories" % len(categories)
print

# split a training set and a test set
y_train, y_test = data_train.target, data_test.target

if opts.labeled_fraction is None:
    fraction = .1
else:
    fraction = opts.labeled_fraction
    if fraction <= 0. or fraction > 1.:
        print "Invalid fraction %.2f"
        sys.exit(1)

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

unlabeled, labeled = split_indices(y_train, fraction)
print "Removing labels of %d random training documents" % len(unlabeled)
print
X_labeled = X_train[labeled]
y_labeled = y_train[labeled]
y_train[unlabeled] = -1

vocabulary = np.array([t for t, i in sorted(vectorizer.vocabulary.iteritems(),
                                            key=itemgetter(1))])


###############################################################################
# Benchmark classifiers
def benchmark(clf, supervised=False):
    print 80 * '_'
    print "Training: "
    print clf
    t0 = time()
    if supervised:
        clf.fit(X_labeled, y_labeled)
    else:
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

        if opts.print_top10:
            print "top 10 keywords per class:"
            for i, category in enumerate(categories):
                top10 = np.argsort(clf.coef_[i, :])[-10:]
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

print 80 * '='
print "Baseline: fully supervised Naive Bayes"
benchmark(MultinomialNB(alpha=.01), supervised=True)
benchmark(BernoulliNB(alpha=.01), supervised=True)

print 80 * '='
print "Naive Bayes trained with Expectation Maximization"
benchmark(SemisupervisedNB(MultinomialNB(alpha=.01)))
benchmark(SemisupervisedNB(BernoulliNB(alpha=.01)))
