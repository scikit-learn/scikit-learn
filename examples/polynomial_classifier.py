"""Polynomial Classifier - An example on sample datasets.

Author: Christoph Hermes <hermes(at)hausmilbe(dot)net>
"""

from sklearn import datasets
from sklearn.polynomial_classifier import PC
from sklearn.metrics import confusion_matrix
from sklearn.decomposition import PCA

import random

def divide_dataset(X,y):
    """Randomly divide dataset given by X,y into train (2/3) and test set (1/3)

    Returns
    -------
    X_test, y_test, X_train, y_train : array-like, shape = [n_samples, n_features]
    """

    n_samples, n_features = X.shape
    p = range(n_samples)
    random.seed(0)
    random.shuffle(p)
    X, y = X[p], y[p]
    third = int(n_samples/3)

    X_test, y_test = X[:third], y[:third]
    X_train, y_train = X[third:], y[third:]

    return X_test, y_test, X_train, y_train

def pc_apply_datasets():
    """Apply Polynomial Classifier to 'iris' and 'digits' datasets and prints
    results on console
    """

    DBs = [('iris', datasets.load_iris()), ('digits', datasets.load_digits())]

    for DB_desc, database in DBs:

        X = database.data
        y = database.target

        X_test, y_test, X_train, y_train = divide_dataset(X,y)

        # PCA: reduce feature dimension by keeping specified amount of
        # information
        pca = PCA(n_components=0.98).fit(X_train)

        # apply classifier
        clf = PC(degree=2)
        y_pred = clf.fit(pca.transform(X_train), y_train).predict(pca.transform(X_test))

        # Compute confusion matrix
        cm = confusion_matrix(y_test, y_pred)

        # print info on console
        print "-" * 80
        print pca, clf
        print "Dataset size: n_samples(Train) = %d, n_samples(Test) = %d" % (X_train.shape[0], X_test.shape[0])
        print "Confusion matrix for '%s' dataset:\n" % DB_desc, cm

if __name__ == '__main__':
    print __doc__
    pc_apply_datasets()

