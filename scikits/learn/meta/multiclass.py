
"""
Multiclass algorithms
======================

This module implements multiclass learning algorithms:
    - one-vs-the-rest
    - one-vs-one
    - error correcting output codes

The algorithms can be used to turn a binary classifier into a multiclass
classifier or to (possibly) improve the accuracy or runtime performance of
multiclass classifiers.
"""

# Author: Mathieu Blondel <mathieu@mblondel.org>
#
# License: BSD Style.

import numpy as np

from scikits.learn.base import BaseEstimator, ClassifierMixin, clone
from scikits.learn.preprocessing import LabelBinarizer
from scikits.learn.metrics.pairwise import euclidean_distances


def fit_binary(estimator, X, y):
    estimator = clone(estimator)
    estimator.fit(X, y)
    return estimator

def predict_binary(estimator, X):
    if hasattr(estimator, "decision_function"):
        return np.ravel(estimator.decision_function(X))
    else:
        # probabilities of the positive class
        return estimator.predict_proba(X)[:, 1]

def fit_ovr(estimator, X, y):
    if not hasattr(estimator, "decision_function") and \
       not hasattr(estimator, "predict_proba"):
        raise ValueError("The base estimator should implement "
                         "decision_function or predict_proba!")

    lb = LabelBinarizer()
    Y = lb.fit_transform(y)
    estimators = [fit_binary(estimator, X, Y[:, i]) for i in range(Y.shape[1])]
    return estimators, lb

def predict_ovr(estimators, label_binarizer, X):
    Y = np.array([predict_binary(e, X) for e in estimators]).T
    return label_binarizer.inverse_transform(Y)


class OneVsRestClassifier(BaseEstimator, ClassifierMixin):

    def __init__(self, estimator):
        self.estimator = estimator

    def fit(self, X, y):
        self.estimators_, self.label_binarizer_ = fit_ovr(self.estimator, X, y)
        return self

    def predict(self, X):
        if not hasattr(self, "estimators_"):
            raise ValueError("The object hasn't been fitted yet!")

        return predict_ovr(self.estimators_, self.label_binarizer_, X)

def fit_ovo_binary(estimator, X, y, i, j):
    cond = np.logical_or(y == i, y == j)
    y = y[cond].copy()
    y[y == i] = 0
    y[y == j] = 1
    ind = np.arange(X.shape[0])
    return fit_binary(estimator, X[ind[cond]], y)

def fit_ovo(estimator, X, y):
    if not hasattr(estimator, "decision_function") and \
       not hasattr(estimator, "predict_proba"):
        raise ValueError("The base estimator should implement "
                         "decision_function or predict_proba!")

    classes = np.unique(y)
    n_classes = classes.shape[0]
    estimators = [fit_ovo_binary(estimator, X, y, classes[i], classes[j])
                    for i in range(n_classes) for j in range(i + 1, n_classes)]

    return estimators, classes

def predict_ovo(estimators, classes, X):
    n_samples = X.shape[0]
    n_classes = classes.shape[0]
    votes = np.zeros((n_samples, n_classes))

    k = 0
    for i in range(n_classes):
        for j in range(i + 1, n_classes):
            pred = estimators[k].predict(X)
            votes[pred == 0, i] += 1
            votes[pred == 1, j] += 1
            k += 1

    return classes[votes.argmax(axis=1)]

class OneVsOneClassifier(BaseEstimator, ClassifierMixin):

    def __init__(self, estimator):
        self.estimator = estimator

    def fit(self, X, y):
        self.estimators_, self.classes_ = fit_ovo(self.estimator, X, y)
        return self

    def predict(self, X):
        if not hasattr(self, "estimators_"):
            raise ValueError("The object hasn't been fitted yet!")

        return predict_ovo(self.estimators_, self.classes_, X)
