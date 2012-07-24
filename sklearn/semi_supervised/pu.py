"""
Algorithms for learning from positive and unlabeled samples
"""

# Author: Lars Buitinck <L.J.Buitinck@uva.nl>
# Copyright 2011-2012, University of Amsterdam
# License: three-clause BSD

from ..base import BaseEstimator, TransformerMixin
from ..utils import atleast2d_or_csr

import numpy as np


def _densify(X):
    return X.toarray() if hasattr(X, 'toarray') else X


class OneDNFTransformer(BaseEstimator, TransformerMixin):
    """1-DNF algorithm for PU learning

    The 1-DNF algorithm is a simple bootstrap method for PU learning on
    multinomially distributed data (e.g., text classification). When fit, it
    estimates a subset of the features to be regarded indicative of positive
    examples (here called the "positive support" in analogy to feature
    selection); when applied, all samples that cannot meet a pre-set threshold
    for each of the positive support features are considered negatives.

    Parameters
    ----------
    thresh : numeric, optional
        Threshold value for transform: samples with all support feature values
        below this number are considered negatives.

    References
    ----------
    B. Liu (2007). Web Data Mining. Springer, pp. 172-173.
    H. Yu, J. Han and K. Chang (2002). PEBL: Positive Example Based Learning
        for web page classification using SVM. Proc. KDD'02, pp. 239-248.
    """
    def __init__(self, thresh=0):
        self.thresh = thresh

    def fit(self, X, y):
        """Estimate which of the unlabeled samples can be considered negatives

        Parameters
        ----------
        X : {array-like, sparse matrix}, shape = [n_samples, n_features]
            Array of samples.

        y : array-like, shape = n_samples
            Label vector; 1 indicates a positive sample, -1 an unlabeled one.
        """
        X = atleast2d_or_csr(X)
        y = np.asanyarray(y)

        # np.ravel to get rid of np.matrix objects
        feature_freq_p = np.ravel(X[np.where(y != -1)[0]].sum(axis=0))
        feature_freq_u = np.ravel(X[np.where(y == -1)[0]].sum(axis=0))

        n_unlabeled = (y == -1).sum()
        feature_freq_u /= n_unlabeled.astype(np.float)
        self.pos_support_ = np.where(feature_freq_p > feature_freq_u)[0]

        return self

    def fit_transform(self, X, y):
        return self.fit(X, y).transform(X, y)

    def transform(self, X, y):
        """Transform {positive, unlabeled} samples to include negatives

        Note that this method, even when handed a sparse matrix, will build a
        temporary dense array of size len(self.pos_support_) * n_unlabeled.

        Parameters
        ----------
        X : {array-like, sparse matrix}, shape = [n_samples, n_features]
            Array of samples.

        y : array-like, shape = n_samples
            Label vector; 1 indicates a positive sample, -1 an unlabeled one.
        """
        X = atleast2d_or_csr(X)
        y = np.array(y)

        # The following computes:
        # for i in the unlabeled X do:
        #   if for any j in pos_support_ s.t. X[i, j] > thresh:
        #       y[i] = 0
        # ... but without loops.

        X_sup = X[np.where(y == -1)[0]][:, self.pos_support_]
        negative = np.all(_densify(X_sup) <= self.thresh, axis=1)

        y[np.where(y == -1)[0][np.where(negative)]] = 0

        return y
