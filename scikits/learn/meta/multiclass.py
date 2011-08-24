
"""
Multiclass algorithms
======================

This module implements multiclass learning algorithms:
    - one-vs-the-rest
    - one-vs-one
    - error correcting output codes
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
    if hasattr(estimators[0], "decision_function"):
        method = "decision_function"
    else:
        method = "predict_proba"

    Y = np.array([getattr(e, method)(X) for e in estimators]).T

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
