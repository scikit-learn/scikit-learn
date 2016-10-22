# Author: Henry Lin <hlin117@gmail.com>

from __future__ import division

import numpy as np
from six.moves import xrange

from sklearn.base import BaseEstimator, TransformerMixin
from sklearn.utils.validation import check_array, check_is_fitted, column_or_1d

class KBinsDiscretizer(BaseEstimator, TransformerMixin):
    """Bins continuous data into k equal width intervals.

    Parameters
    ----------
    n_bins : int (default=2)
        The number of bins to produce. The intervals for the bins are
        determined by the minimum and maximum of the input data.
        Raises ValueError if n_bins < 2.

    categorical_features : int array-like (default=None)
        Column indexes of categorical features. Categorical features
        will not be discretized, and will remain as is after the transform.

    Attributes
    ----------
    mins_ : float array, shape (n_features_,)
        Minimum value per feature in the input data.

    maxes_ : float array, shape (n_features,)
        Maximum value per feature in the input data.

    n_features_ : int
        Number of features in X.

    cut_points_ : int array, shape (n_bins - 1, n_features_)
        Each column of cut_points_ is a set of cut points which describe
        how the values of the corresponding feature are binned. For a pair
        of adjacent values in a column, a, b, values between [a, b) will
        have the same discretization value.

    Example
    -------
    >>> X = [[-2, 1, -4,   -1], \
             [-1, 2, -3, -0.5], \
             [ 0, 3, -2,  0.5], \
             [ 1, 4, -1,    2]]
    >>> dis = KBinsDiscretizer(n_bins=3)
    >>> dis.fit(X) # doctest: +ELLIPSIS, +NORMALIZE_WHITESPACE
    KBinsDiscretizer(...)
    >>> dis.cut_points_ # doctest: +NORMALIZE_WHITESPACE
    array([[-1., 2., -3., 0.],
           [ 0., 3., -2., 1.]])
    >>> dis.transform(X) # doctest: +ELLIPSIS, +NORMALIZE_WHITESPACE
    array([[0, 0, 0, 0],
           [1, 1, 1, 0],
           [2, 2, 2, 1],
           [2, 2, 2, 2]])

    Categorical features can be specified in the constructor.

    >>> dis = KBinsDiscretizer(n_bins=3, categorical_features=[1])
    >>> dis.fit(X) # doctest: +ELLIPSIS, +NORMALIZE_WHITESPACE
    KBinsDiscretizer(...)

    For each feature mentioned as a categorical feature, this feature
    will be represented as nan values in the discretizer's cut points.
    >>> dis.cut_points_ # doctest: +NORMALIZE_WHITESPACE
    array([[ -1., nan, -3., 0.],
           [  0., nan, -2., 1.]])

    The discretization process will not touch any column marked as categorical.
    >>> dis.transform(X) # doctest: +ELLIPSIS, +NORMALIZE_WHITESPACE
    array([[ 0., 1., 0., 0.],
           [ 1., 2., 1., 0.],
           [ 2., 3., 2., 1.],
           [ 2., 4., 2., 2.]])
    """


    def __init__(self, n_bins=2, categorical_features=None):
        self.n_bins = n_bins
        self.categorical_features = categorical_features

        # Attributes
        self.mins_ = None
        self.maxes_ = None
        self.n_features_ = None
        self.cut_points_ = None
        self.categorical_features_set_ = None


    def fit(self, X, y=None):
        """Fits the estimator.

        Parameters
        ----------
        X : numeric array-like, shape (n_samples, n_features)
            The input data containing features to be discretized. Indexes
            of categorical features are indicated by categorical_features
            from the constructor.

        Returns
        -------
        self
        """

        if self.n_bins < 2:
            raise ValueError("Discretizer received an invalid number "
                             "of bins. Received {0}, expected at least 2."
                             .format(self.n_bins))

        # TODO: Determine behavior for 1D input
        X = check_array(X, dtype=float)

        self._check_categorical_features(self.categorical_features, X.shape[1])

        self._fit(X)
        return self

    @staticmethod
    def _check_categorical_features(cat_features, n_features):
        """Only allow integer indices, or None.
        """
        if cat_features is None:
            return

        cat_features = check_array(cat_features, ensure_2d=False, dtype=int)
        cat_features = column_or_1d(cat_features)

        if len(set(cat_features)) != cat_features.shape[0]:
            raise ValueError("Duplicate categorical column indexes found.")

        if np.all(cat_features >= 0) and np.all(cat_features < n_features):
            return

        raise ValueError("Invalid categorical feature index.")

    def _cut_points_from_index(self, col_idx):
        """Determines the cut points of a feature based on the column index.
        """
        n_bins = self.n_bins
        if col_idx in self.categorical_features_set_:
            return [np.nan] * (self.n_bins - 1)

        _min = self.mins_[col_idx]
        _max = self.maxes_[col_idx]
        spacing = (_max - _min) / n_bins
        return np.linspace(_min + spacing, _max - spacing, num=n_bins - 1)


    def _fit(self, X):
        """Fits the attributes of the discretizer.
        """
        n_features = X.shape[1]
        if self.categorical_features:
            self.categorical_features_set_ = set(self.categorical_features)
        else:
            self.categorical_features_set_ = {}

        self.mins_ = np.min(X, axis=0)
        self.maxes_ = np.max(X, axis=0)

        cut_points = [self._cut_points_from_index(i)
                      for i in xrange(n_features)]

        self.cut_points_ =  np.array(cut_points).T
        self.n_features_ = n_features


    def transform(self, X, y=None):
        """Discretizes the data. Leaves any categorical features as is.

        Parameters
        ----------
        X : numeric array-like, shape (n_samples, n_features)
            Data containing features to be discretized.

        Returns
        -------
        T : numeric array-like, shape (n_samples, n_features)
            Categorical features will not be transformed. If all features
            are continuous, T will be an int array. Otherwise T will be
            a float array.

        """
        check_is_fitted(self, ["mins_", "maxes_", "n_features_", "cut_points_"])
        X = check_array(X, dtype=float)

        if X.shape[1] != self.n_features_:
            raise ValueError("Transformed array does not have correct number "
                             "of features. Expecting {0}, received {1}."
                             .format(self.n_features_, X.shape[1]))

        return self._transform(X)

    def _transform_from_index(self, X, col_idx):
        X_col = X[:, col_idx]
        if col_idx in self.categorical_features_set_:
            return X_col
        return np.digitize(X_col, self.cut_points_[:, col_idx])


    def _transform(self, X):
        """Discretizes the continuous features into k integers.
        """
        discretized = [self._transform_from_index(X, col_idx)
                       for col_idx in xrange(self.n_features_)]

        return np.array(discretized).T


    def _transform_deprecated(self, X):
        """Discretizes the continuous features into k integers.

        NOTE: Not used until further notice. Under consideration for sparse case.

        The process used below preserves the sparsity of the data by
        transforming the data, and rounding the data down to the
        nearest integer. These integers represent the bins sought.

        The process is as follows.
        1. Let X_trans = X * n_bins / (X.max() - X.min()).
            Observe that X_trans is a dataset of width k.
        2. epsilon = abs(X_trans.min() - ceil(X_trans.min())).
            episilon < 1 is the distance of X.min() to ceil(X.min())
        3. X_trans[X_trans != 0] += epsilon
            Shifts the data. X_trans.min() is now an integer. Note this
            operation can be done on sparse matrices.
        4. output = floor(X_trans)
            The output now contains integers which represent the
            discretized values. Observe that even if we added epsilon to
            zero values in step 3, these zeros will still be discretized
            to zero.
        """
        # Rescale data
        X_trans = X / (self.max_ - self.min_)
        X_trans *= self.n_bins

        # Shift data
        epsilon = np.abs(X_trans.min() - np.ceil(X_trans.min()))
        X_trans[X_trans != 0] += epsilon

        # Generate output
        discretized = np.floor(X_trans)

        # Corner case arises when maximum values across each column of X
        # get shifted to the next integer.
        # TODO: Use numpy clip here
        for col in discretized.T:
            col_min = col.min()
            col_max = col.max()
            if col_max - col_min > self.n_bins - 1:
                col[col == col_max] = col_max - 1

        return discretized
