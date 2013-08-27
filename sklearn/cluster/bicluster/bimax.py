"""Implements the BiMax biclustering algorithm.

Authors : Kemal Eren
License: BSD 3 clause

"""
from abc import ABCMeta

import numpy as np

from sklearn.base import BaseEstimator, BiclusterMixin
from sklearn.externals import six

from .utils import get_indicators

from ._biclique import find_bicliques


class BiMax(six.with_metaclass(ABCMeta, BaseEstimator,
                               BiclusterMixin)):
    """Method to find all maximal biclusters in a boolean array."""

    def __init__(self):
        pass

    def fit(self, X):
        """Creates a biclustering for X.

        Parameters
        ----------
        X : array-like, shape (n_samples, n_features)

        """
        # TODO: check X
        result = list(find_bicliques(X))

        all_rows = []
        all_cols = []
        n_rows = X.shape[0]
        for nodes in result:
            rows = list(n for n in nodes if n < n_rows)
            cols = list(n - n_rows for n in nodes if n >= n_rows)
            if not rows or not cols:
                continue
            rows, cols = np.array(rows), np.array(cols)
            row_idx, col_idx = get_indicators(rows, cols, X.shape)
            all_rows.append(row_idx)
            all_cols.append(col_idx)
        self.rows_ = np.vstack(all_rows)
        self.columns_ = np.vstack(all_cols)
