# -*- coding: utf-8 -*-
"""χ² (chi-square) feature selection"""

import numpy as np
from scipy.stats import chisquare

from ..preprocessing import LabelBinarizer
from ..utils import safe_asanyarray
from ..utils.extmath import safe_sparse_dot


class Chi2(BaseEstimator, TransformerMixin):
    """Select best features by the χ² statistic

    This transformer can be used to select the n_features features with the
    highest values for the χ² (chi-square) statistic from multinomially
    distributed data (e.g., term counts in document classification).

    Note that this class does not perform a significance test.

    Attributes
    ----------
    top_features_ : array, dtype = int, shape = [n_features]
        The top n_features features by χ² value.

    References
    ----------
    Y. Yang and J.P. Pedersen (1997). A comparative study on feature selection
        in text categorization. Proc. ICML, pp. 412-420.
        http://nyc.lti.cs.cmu.edu/yiming/Publications/icml97.ps
    """

    def __init__(self, n_features):
        self.n_features = n_features

    def fit(self, X, y):
        """Find the best features in a training set."""

        # XXX: we might want to do some of the following in logspace instead
        # for numerical stability.
        X = safe_asanyarray(X)
        Y = LabelBinarizer().fit_transform(y)
        if Y.shape[1] == 1:
            Y = np.concatenate((1 - Y, Y), axis=1)

        observed = safe_sparse_dot(Y.T, X)      # n_classes * n_features

        #N = observed.sum(dtype=np.float64)
        feature_count = np.atleast_2d(X.sum(axis=0))
        class_prob = np.atleast_2d(Y.sum(axis=0) / Y.sum())
        expected = safe_sparse_dot(class_prob.T, feature_count)

        # We lied in the docstring: we do perform a significance test,
        # only to then throw the result away.
        ch2, _ = chisquare(observed, expected)

        self.top_features_ = np.argsort(ch2)[-self.n_features:]
        return self

    def transform(self, X, y=None):
        """Pass through the features of X with the highest χ² in training."""
        X = safe_asanyarray(X)
        return X[:, self.top_features_]
