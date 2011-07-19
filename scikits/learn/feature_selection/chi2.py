# -*- coding: utf-8 -*-
"""χ² (chi-square) feature selection"""

import numpy as np
from scipy.stats import chisquare

from ..base import BaseEstimator, TransformerMixin
from ..preprocessing import LabelBinarizer
from ..utils import safe_asanyarray
from ..utils.extmath import safe_sparse_dot


class Chi2(BaseEstimator, TransformerMixin):
    """Select best features by the χ² statistic

    This transformer can be used to select the n_features features with the
    highest values for the χ² (chi-square) statistic from multinomially
    distributed data (e.g., term counts in document classification) relative
    to the classes.

    Recall that the χ² statistic measures dependence between stochastic
    variables, so this transformer "weeds out" the features that are the most
    likely to be independent of class and therefore irrelevant for
    classification.

    Parameters
    ----------
    n_features : int
        Number of features to select.

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
        """Find the best features in a training set.

        Parameters
        ----------
        X : {array-like, sparse matrix}, shape = [n_samples, n_features_in]
            Sample vectors.

        y : array-like, shape = n_samples
            Target vector (class labels).

        Returns
        -------
        self
        """

        # XXX: we might want to do some of the following in logspace instead
        # for numerical stability.
        X = safe_asanyarray(X)
        Y = LabelBinarizer().fit_transform(y)
        if Y.shape[1] == 1:
            Y = np.concatenate((1 - Y, Y), axis=1)

        observed = safe_sparse_dot(Y.T, X)      # n_classes * n_features

        feature_count = np.atleast_2d(X.sum(axis=0))
        feature_count = np.asarray(feature_count)   # stupid numpy.matrix!
        class_prob = np.atleast_2d(Y.sum(axis=0) / Y.sum())
        expected = safe_sparse_dot(class_prob.T, feature_count)

        # The p-values (probabilities of independence) are monotonically
        # decreasing in χ², so we need only one of both values.
        ch2, _ = chisquare(observed, expected)

        self.top_features_ = np.argsort(ch2)[-self.n_features:]
        return self

    def transform(self, X, y=None):
        """Pass through the features of X with the highest χ² in training.

        Parameters
        ----------
        X : {array-like, sparse matrix}, shape = [n_samples, n_features_in]
            Sample vectors. If X is a sparse matrix, it must be in a format
            that supports "fancy indexing" (slicing).

        Returns
        -------
        Xsel : {array-like, sparse matrix}, shape = [n_samples, n_features]
            Sample vectors with all but the top n_features columns removed.
        """
        X = safe_asanyarray(X)
        return X[:, self.top_features_]
