"""
======================================================
Classification of text documents using sparse features
======================================================

This is an example showing how scikit-learn can be used to classify documents
by topics using a bag-of-words approach. This example uses a scipy.sparse
matrix to store the features and demonstrates various classifiers that can
efficiently handle sparse matrices.

The dataset used in this example is the 20 newsgroups dataset. It will be
automatically downloaded, then cached.

The bar plot indicates the accuracy, training time (normalized) and test time
(normalized) of each classifier.

"""

# Author: Peter Prettenhofer <peter.prettenhofer@gmail.com>
#         Olivier Grisel <olivier.grisel@ensta.org>
#         Mathieu Blondel <mathieu@mblondel.org>
#         Lars Buitinck
# License: BSD 3 clause

from __future__ import print_function
from __future__ import division

import numpy as np
import matplotlib.pyplot as plt

import logging
import sys
from time import time
from optparse import OptionParser

from sklearn.datasets import fetch_20newsgroups
from sklearn.ensemble import RandomForestClassifier
from sklearn.feature_extraction.text import TfidfVectorizer, HashingVectorizer
from sklearn.feature_selection import SelectFromModel, SelectKBest, chi2
from sklearn.svm import LinearSVC
from sklearn.linear_model import (SGDClassifier, Perceptron, RidgeClassifier,
                                  PassiveAggressiveClassifier)
from sklearn.naive_bayes import BernoulliNB, MultinomialNB
from sklearn.neighbors import KNeighborsClassifier, NearestCentroid
from sklearn.pipeline import Pipeline
from sklearn.utils.extmath import density
from sklearn import metrics


# Display progress logs on stdout
logging.basicConfig(level=logging.INFO,
                    format='%(asctime)s %(levelname)s %(message)s')


# Parse command line arguments
op = OptionParser()
op.add_option('--all_categories',
              action='store_true', dest='all_categories',
              help='Whether to use all categories or not.')
op.add_option('--confusion_matrix',
              action='store_true', dest='print_cm',
              help='Print the confusion matrix.')
op.add_option('--filtered',
              action='store_true',
              help='Remove newsgroup information that is easily overfit: '
                   'headers, signatures, and quoting.')
op.add_option('--report',
              action='store_true', dest='print_report',
              help='Print a detailed classification report.')
op.add_option('--top10',
              action='store_true', dest='print_top10',
              help='Print ten most discriminative terms per class '
                   'for every classifier.')
op.add_option('--use_hashing',
              action='store_true',
              help='Use a hashing vectorizer.')
op.add_option('--chi2_select',
              action='store', type=int, dest='select_chi2',
              help='Select some number of features using a '
                   'chi-squared test.')
op.add_option('--n_features',
              action='store', type=int, default=2 ** 16,
              help='n_features when using the hashing vectorizer.')


def benchmark(clf):
    """Trains and evaluates a classifier.
    Args:
        clf: classifier instance
    Returns:
        clf description, accuracy score, train & test times
    """
    print('_' * 80)
    print('Training: ')
    print(clf)

    t0 = time()
    clf.fit(X_train, y_train)
    train_time = time() - t0
    print('Train Time: %0.3fs' % train_time)

    t0 = time()
    pred = clf.predict(X_test)
    test_time = time() - t0
    print('Test Time:  %0.3fs' % test_time)

    score = metrics.accuracy_score(y_test, pred)
    print('Accuracy:   %0.3f' % score)

    if hasattr(clf, 'coef_'):
        print('Dimensionality: %d' % clf.coef_.shape[1])
        print('Density:        %f' % density(clf.coef_))

        if opts.print_top10 and feature_names is not None:
            print('Top 10 Keywords per Class:')
            for i, label in enumerate(target_names):
                top10 = np.argsort(clf.coef_[i])[-10:]
                print(trim('%s: %s' % (label,
                                       ' '.join(feature_names[top10]))))
        print()

    if opts.print_report:
        print('Classification Report:')
        print(metrics.classification_report(y_test, pred,
                                            target_names=target_names))

    if opts.print_cm:
        print('Confusion Matrix:')
        print(metrics.confusion_matrix(y_test, pred))

    print()
    clf_descr = str(clf).split('(')[0]

    return clf_descr, score, train_time, test_time


def size_mb(docs):
    """Returns size of all documents in MB"""
    return sum(len(s.encode('utf-8')) for s in docs) / 1e6


def trim(s):
    """Trim string to fit on terminal (assuming 80-column display)"""
    return s if len(s) <= 80 else s[:77] + '...'


def is_interactive():
    return not hasattr(sys.modules['__main__'], '__file__')

# work-around for Jupyter notebook and IPython console
argv = [] if is_interactive() else sys.argv[1:]
(opts, args) = op.parse_args(argv)
if len(args) > 0:
    op.error('Received invalid arguments: %s' % args)
    op.print_help()
    sys.exit(1)

print(__doc__)
op.print_help()
print()


###############################################################################
# Load some categories from the training set
if opts.all_categories:
    categories = None
else:
    categories = [
        'alt.atheism',
        'talk.religion.misc',
        'comp.graphics',
        'sci.space',
    ]

if opts.filtered:
    remove = ('headers', 'footers', 'quotes')
else:
    remove = ()

print('Loading 20 newsgroups dataset for categories:')
print(categories if categories else 'all')

data_train = fetch_20newsgroups(subset='train', categories=categories,
                                shuffle=True, random_state=42,
                                remove=remove)

data_test = fetch_20newsgroups(subset='test', categories=categories,
                               shuffle=True, random_state=42,
                               remove=remove)
print('Data loaded.')

# Order of labels in `target_names` can be different from `categories`
target_names = data_train.target_names

data_train_size_mb = size_mb(data_train.data)
data_test_size_mb = size_mb(data_test.data)

print('%d documents - %0.3fMB (training set)' % (
    len(data_train.data), data_train_size_mb))
print('%d documents - %0.3fMB (test set)' % (
    len(data_test.data), data_test_size_mb))
print('%d categories' % len(categories))
print()

# Split a training set and a test set
y_train, y_test = data_train.target, data_test.target

print('Extracting features from the training data using a '
      'sparse %s vectorizer.' % ('hashing' if opts.use_hashing else 'tf-idf'))
t0 = time()
if opts.use_hashing:
    vectorizer = HashingVectorizer(stop_words='english', non_negative=True,
                                   n_features=opts.n_features)
    X_train = vectorizer.transform(data_train.data)
else:
    vectorizer = TfidfVectorizer(sublinear_tf=True, max_df=0.5,
                                 stop_words='english')
    X_train = vectorizer.fit_transform(data_train.data)

duration = time() - t0
print('Done in %fs at %0.3fMB/s' % (duration, data_train_size_mb / duration))
print('n_samples: %d, n_features: %d' % X_train.shape)
print()

print('Extracting features from the test data using the same vectorizer')
t0 = time()
X_test = vectorizer.transform(data_test.data)
duration = time() - t0
print('Done in %fs at %0.3fMB/s' % (duration, data_test_size_mb / duration))
print('n_samples: %d, n_features: %d' % X_test.shape)
print()

# Mapping from integer feature name to original token string
if opts.use_hashing:
    feature_names = None
else:
    feature_names = vectorizer.get_feature_names()

if opts.select_chi2:
    print('Extracting %d best features by a chi-squared test.' %
          opts.select_chi2)
    t0 = time()
    ch2 = SelectKBest(chi2, k=opts.select_chi2)
    X_train = ch2.fit_transform(X_train, y_train)
    X_test = ch2.transform(X_test)

    if feature_names:
        # keep selected feature names
        feature_names = [feature_names[i] for i
                         in ch2.get_support(indices=True)]

    print('Done in %fs' % (time() - t0))
    print()

if feature_names:
    feature_names = np.asarray(feature_names)


###############################################################################
# Benchmark classifiers
results = []
for clf, name in (
        (RidgeClassifier(tol=1e-2, solver='sag'), 'Ridge Classifier'),
        (Perceptron(n_iter=50), 'Perceptron'),
        (PassiveAggressiveClassifier(n_iter=50), 'Passive-Aggressive'),
        (KNeighborsClassifier(n_neighbors=10), 'kNN'),
        (RandomForestClassifier(n_estimators=100), 'Random Forest'),
        (SGDClassifier(alpha=.0001, n_iter=50, penalty='elasticnet'),
            'SGD (Elastic-Net Penalty)'),
        (NearestCentroid(), 'NearestCentroid (aka Rocchio classifier)')):
    print('=' * 80)
    print(name)
    results.append(benchmark(clf))

for penalty in ['l2', 'l1']:
    print('=' * 80)
    print('%s Penalty' % penalty.upper())
    print('----------')

    # Liblinear and SGD models
    for clf, name in (
            (LinearSVC(penalty=penalty, dual=False, tol=1e-3), 'LibLinear'),
            (SGDClassifier(alpha=.0001, n_iter=50, penalty=penalty), 'SGD')):
        print('%s (%s)' % (name, penalty.upper()))
        results.append(benchmark(clf))

# Train sparse Naive Bayes classifiers
print('=' * 80)
print('Naive Bayes')
print('-----------')
for clf, name in (
        (MultinomialNB(alpha=.01),  'Multinomial NB'),
        (BernoulliNB(alpha=.01), 'Bernoulli NB')):
    print(name)
    results.append(benchmark(clf))

print('=' * 80)
print('LinearSVC with L1-based feature selection')
# The smaller C, the stronger the regularization.
# The more regularization, the more sparsity.
results.append(benchmark(Pipeline([
  ('feature_selection', SelectFromModel(LinearSVC(penalty='l1', dual=False,
                                                  tol=1e-3))),
  ('classification', LinearSVC(penalty='l2'))])))

# Make some plots
# Compute plot indices for results
indices = np.arange(len(results))

# Sort results by score (accuracy)
results = [[result[i] for result in sorted(results, key=lambda r: r[1])]
           for i in range(4)]

clf_names, score, training_time, test_time = results
training_time = np.array(training_time) / np.max(training_time)
test_time = np.array(test_time) / np.max(test_time)

plt.figure(figsize=(12, 8))
plt.title('Classifier Accuracy')
plt.barh(indices, score, .2, label='score', color='navy')
plt.barh(indices + .3, training_time, .2, label='training time',
         color='c')
plt.barh(indices + .6, test_time, .2, label='test time', color='darkorange')
plt.yticks(())
plt.legend(loc='best')
plt.subplots_adjust(left=.25)
plt.subplots_adjust(top=.95)
plt.subplots_adjust(bottom=.05)

for i, c in zip(indices, clf_names):
    plt.text(-.3, i, c)

plt.show()
