"""Polynomial Classifier - An example on sample datasets.

Author: Christoph Hermes <hermes(at)hausmilbe(dot)net>
"""

from sklearn import datasets
from sklearn.polynomial_classifier import PolynomialClassifier
from sklearn.metrics import confusion_matrix
from sklearn.decomposition import PCA

import random

print __doc__

# list of databases the Polynomial Classifier is applied to
DBs = [('iris', datasets.load_iris()), ('digits', datasets.load_digits())]

for DB_desc, DB_data in DBs:

    X = DB_data.data
    y = DB_data.target

    # Randomly divide dataset into train (2/3) and test set (1/3)
    n_samples, n_features = X.shape
    p = range(n_samples)
    random.seed(0)
    random.shuffle(p)
    X, y = X[p], y[p]
    third = int(n_samples / 3)
    X_test, y_test = X[:third], y[:third]
    X_train, y_train = X[third:], y[third:]

    # PCA: reduce feature dimension by keeping specified amount of information
    pca = PCA(n_components=0.95).fit(X_train)

    # apply classifier
    clf = PolynomialClassifier(degree=2)
    X_train_pca = pca.transform(X_train)
    X_test_pca = pca.transform(X_test)
    y_pred = clf.fit(X_train_pca, y_train).predict(X_test_pca)

    # Compute confusion matrix
    cm = confusion_matrix(y_test, y_pred)

    # print info on console
    print "-" * 80
    print pca, clf
    print "Dataset size: n_samples(Train) = %d, n_samples(Test) = %d" % (
        X_train.shape[0], X_test.shape[0])
    print "Confusion matrix for '%s' dataset:\n" % DB_desc, cm

