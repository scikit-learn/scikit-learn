from __future__ import division

__author__ = "Henry Lin <hlin117@gmail.com>"

import numpy as np
import scipy.sparse as sp
from sklearn.base import BaseEstimator, TransformerMixin
from sklearn.utils import column_or_1d, check_array
from sklearn.utils.validation import check_is_fitted

__all__ = [
    "Discretizer"
]


class Discretizer(BaseEstimator, TransformerMixin):
    """Bins continuous data into k equal width intervals.

    Parameters
    ----------
    n_bins : int (default=2)
        The number of bins to produce. The intervals for the bins are
        determined by the minimum and maximum of the input data.
        Raises ValueError if n_bins < 2.

    categorical_features : array-like or None (default=None)
        Specifies the indices of categorical columns which are not to
        be discretized. If None, assumes that all columns are continuous
        features.

    Attributes
    ----------
    min_ : float
        The minimum value of the input data.

    max_ : float
        The maximum value of the input data.

    cut_points_ : array, shape [numBins - 1]
        Contains the boundaries for which the data lies. Each interval
        has an open left boundary, and a closed right boundary.

    n_features_ : int
        The number of features from the original dataset.

    continuous_mask_ : array, shape [n_continuous_features_]
        An array of booleans, representing columns for which there are
        continuous features.

    Example
    -------
    >>> # X has two examples, with X[:, 2] as categorical
    >>> X = [[0, 1, 0, 4], \
             [6, 7, 1, 5]]
    >>> from sklearn.preprocessing import Discretizer
    >>> discretizer = Discretizer(n_bins=3, categorical_features=[2])
    >>> discretizer.fit(X) # doctest: +ELLIPSIS, +NORMALIZE_WHITESPACE
    Discretizer(...)
    >>> discretizer.cut_points_ # doctest: +ELLIPSIS, +NORMALIZE_WHITESPACE
    array([[ 2. , 3. , 4.3333...],
           [ 4. , 5. , 4.6666...]])
    >>> # Transforming X will move the categorical features to the last indices
    >>> discretizer.transform(X)
    array([[0, 0, 0, 0],
           [2, 2, 2, 1]])
    """
    sparse_formats = ['csr', 'csc']

    def __init__(self, n_bins=2, categorical_features="all"):
        if n_bins < 2:
            raise ValueError("Discretizer received an invalid number "
                             "of bins")
        self.n_bins = n_bins
        self.categorical_features = categorical_features

        # Attributes
        self.min_ = None
        self.max_ = None
        self.cut_points_ = None
        self.n_features_ = None
        self.continuous_mask_ = None

    def _check_sparse(self, X, ravel=True):
        if ravel:
            return X.toarray().ravel() if sp.issparse(X) else X
        return X.toarray() if sp.issparse(X) else X

    def _set_continuous_mask(self):
        """Sets a boolean array that determines which columns are
        continuous features.
        """
        if self.categorical_features == "all":
            self.n_continuous_features_ = self.n_features_
            self.continuous_mask_ = np.ones(self.n_features_, dtype=bool)
            return

        if len(self.categorical_features) > self.n_features_:
            raise ValueError("The number of categorical indices is more than "
                             "the number of features in the dataset.")

        self.continuous_mask_ = np.ones(self.n_features_, dtype=bool)
        self.continuous_mask_[self.categorical_features] = False

        # Checks if there are duplicate columns from self.categorical_features,
        # such as [-1, 0, 0, self.n_features - 1]
        if self.n_features_ - len(self.categorical_features) \
                != self.continuous_mask_.sum():
            raise ValueError("Duplicate indices detected from input "
                             "categorical indices.")

    def fit(self, X, y=None):
        """Finds the intervals of interest from the input data.

        Parameters
        ----------
        X : {array-like, sparse-matrix}, shape [n_samples, n_features]
            The input data containing continuous features to be
            discretized. Categorical columns are indicated by
            categorical_columns from the constructor.

        Returns
        -------
        self
        """

        X = check_array(X, accept_sparse=self.sparse_formats)
        self.n_features_ = X.shape[1]

        # Set the mask for the continuous features in the array
        self._set_continuous_mask()

        continuous = X[:, self.continuous_mask_]

        self.min_ = self._check_sparse(continuous.min(axis=0))
        self.max_ = self._check_sparse(continuous.max(axis=0))
        cut_points = list()
        for min_, max_ in zip(self.min_, self.max_):
            points = np.linspace(min_, max_, num=self.n_bins, endpoint=False)[1:]
            cut_points.append(points.reshape(-1, 1))
        self.cut_points_ = np.hstack(cut_points)
        return self

    def transform(self, X, y=None):
        """Discretizes the input data.

        Parameters
        ----------
        X : {array-like, sparse-matrix}, shape [n_samples, n_features]
            The input data containing continuous features to be
            discretized. Categorical columns are indicated by
            categorical_columns from the constructor.

        Returns
        -------
        T : array-like, shape [n_samples, n_features]
            Array containing discretized features alongside the
            categorical features. Categorical features are represented
            as the last features of the output.

        """
        check_is_fitted(self, ["cut_points_", "continuous_mask_"])
        X = check_array(X, accept_sparse=self.sparse_formats)
        if X.shape[1] != self.n_features_:
            raise ValueError("Transformed array does not have currect number "
                             "of features. Expecting {0}, received {1}"
                             .format(self.n_features_, X.shape[1]))

        continuous = X[:, self.continuous_mask_]
        discretized = list()
        for cut_points, cont in zip(self.cut_points_.T, continuous.T):
            cont = self._check_sparse(cont)  # np.searchsorted can't handle sparse
            dis_features = np.searchsorted(cut_points, cont)
            discretized.append(dis_features.reshape(-1, 1))

        discretized = np.hstack(discretized)
        categorical = self._check_sparse(X[:, ~self.continuous_mask_], ravel=False)
        return np.hstack((discretized, categorical))

