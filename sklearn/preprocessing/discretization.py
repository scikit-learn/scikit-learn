# Author: Henry Lin <hlin117@gmail.com>

from __future__ import division, absolute_import

import numbers
import numpy as np

from sklearn.base import BaseEstimator, TransformerMixin
from sklearn.utils.validation import check_array, check_is_fitted, column_or_1d


class KBinsDiscretizer(BaseEstimator, TransformerMixin):
    """Bins continuous data into k equal width intervals.

    Parameters
    ----------
    n_bins : int or array-like, shape (n_features_,) (default=2)
        The number of bins to produce. The intervals for the bins are
        determined by the minimum and maximum of the input data.
        Raises ValueError if n_bins < 2.

        If n_bins is an array, and there is a categorical feature at index i,
        n_bins[i] will be ignored.

    categorical_features : int array-like (default=None)
        Column indices of categorical features. Categorical features
        will not be discretized and will be ignored by the transformation.

    Attributes
    ----------
    categorical_features_ : int array, shape (n_features_,)
        Indices of categorical features. These are not affected by the
        transformation.

    min_ : float array, shape (n_features_,)
        Minimum value per feature in X. A categorical feature at index i will
        have min_[i] == 0.

    ptp_ : float array, shape (n_features_,)
        X.max(axis=0) - X.min(axis=0). A categorical feature at index i will
        have ptp_[i] == 1.

    clip_min_ : float array, shape (n_features_,)
        Ensures that continuous features achieve a minimum value of 0 after
        discretization. A categorical feature at index i will have
        clip_min_[i] == -np.inf.

    clip_max_ : float array, shape (n_features_,)
        Ensures that continuous features achieve a maximum value of n_bins - 1
        after discretization. A categorical feature at index i will have
        clip_max_[i] == np.inf.

    n_bins_ : int array, shape (n_features_,)
        Number of bins per feature. A categorical feature at index i will
        have n_bins_[i] == 1.

    n_features_ : int
        Number of features in X.

    Example
    -------
    >>> X = [[-2, 1, -4,   -1], \
             [-1, 2, -3, -0.5], \
             [ 0, 3, -2,  0.5], \
             [ 1, 4, -1,    2]]
    >>> dis = KBinsDiscretizer(n_bins=3)
    >>> dis.fit(X) # doctest: +ELLIPSIS, +NORMALIZE_WHITESPACE
    KBinsDiscretizer(...)
    >>> dis.transform(X) # doctest: +NORMALIZE_WHITESPACE
    array([[ 0., 0., 0., 0.],
           [ 1., 1., 1., 0.],
           [ 2., 2., 2., 1.],
           [ 2., 2., 2., 2.]])

    Categorical features can be specified in the constructor.

    >>> dis = KBinsDiscretizer(n_bins=3, categorical_features=[1])
    >>> dis.fit(X) # doctest: +ELLIPSIS, +NORMALIZE_WHITESPACE
    KBinsDiscretizer(...)

    Discretization will ignore categorical features.
    >>> dis.transform(X) # doctest: +NORMALIZE_WHITESPACE
    array([[ 0., 1., 0., 0.],
           [ 1., 2., 1., 0.],
           [ 2., 3., 2., 1.],
           [ 2., 4., 2., 2.]])
    """

    def __init__(self, n_bins=2, categorical_features=None):
        self.n_bins = n_bins
        self.categorical_features = categorical_features

    def fit(self, X, y=None):
        """Fits the estimator.

        Parameters
        ----------
        X : numeric array-like, shape (n_samples, n_features)
            Data to be discretized.

        Returns
        -------
        self
        """
        X = check_array(X, dtype=float, ensure_2d=False)

        if X.ndim == 1:
            raise ValueError("Reshape your data.")

        n_features = X.shape[1]
        self._check_categorical_features(self.categorical_features, X.shape[1])
        self._check_n_bins(self.n_bins, n_features)
        self.n_features_ = n_features

        if self.categorical_features:
            self.categorical_features_ = np.asarray(self.categorical_features)
        else:
            self.categorical_features_ = np.array([], dtype='int64')

        # Manipulating min and max so categorical data is not affected by
        # transformation.
        min = np.min(X, axis=0)
        max = np.max(X, axis=0)

        min[self.categorical_features_] = 0
        max[self.categorical_features_] = 1

        self.min_ = min
        self.ptp_ = max - min

        # Set n_bins_ so categorical data can be ignored by transformation.
        if isinstance(self.n_bins, numbers.Number):
            n_bins = np.ones(self.n_features_) * self.n_bins
        else:
            n_bins = check_array(self.n_bins, dtype=int, ensure_2d=False,
                                 copy=True)

        n_bins[self.categorical_features_] = 1
        self.n_bins_ = n_bins

        # Clipping features ensures outliers are within boundary after
        # transformation, but also ensures categorical features can be ignored.
        numeric_features = np.delete(np.arange(self.n_features_),
                                     self.categorical_features_)
        clip_min = np.repeat(-np.inf, self.n_features_)
        clip_min[numeric_features] = 0

        clip_max = np.repeat(np.inf, self.n_features_)
        clip_max[numeric_features] = self.n_bins_[numeric_features] - 1

        self.clip_min_ = clip_min
        self.clip_max_ = clip_max

        same_min_max = np.where(self.ptp_ == 0)[0]
        if len(same_min_max) > 0:
            import warnings
            warnings.warn("Fitted X contained continuous features at indices "
                          "{}. All of these features will be discretized to "
                          "the value 0."
                          .format(", ".join(str(i) for i in same_min_max)))

        return self

    @staticmethod
    def _check_n_bins(n_bins, n_features):
        if isinstance(n_bins, numbers.Number):
            if n_bins < 2:
                raise ValueError("Discretizer received an invalid number "
                                 "of bins. Received {0}, expected at least 2."
                                 .format(n_bins))
            return

        n_bins = check_array(n_bins, dtype=int, ensure_2d=False)

        if n_bins.ndim > 1 or n_bins.shape[0] != n_features:
            raise ValueError("n_bins must be a scalar or array "
                             "of shape (n_features_,).")

        violating_indices = np.arange(n_features)[n_bins < 2]
        if len(violating_indices) > 0:
            indices = ", ".join(str(i) for i in violating_indices)
            raise ValueError("Discretizer received an invalid number "
                             "of bins at indices {}. Number of bins "
                             "must be at least 2."
                             .format(indices))

    @staticmethod
    def _check_categorical_features(cat_features, n_features):
        if cat_features is None:
            return

        cat_features = check_array(cat_features, ensure_2d=False, dtype=int)
        cat_features = column_or_1d(cat_features)

        if len(set(cat_features)) != cat_features.shape[0]:
            raise ValueError("Duplicate categorical column indices found.")

        if np.all(cat_features >= 0) and np.all(cat_features < n_features):
            return

        raise ValueError("Invalid categorical feature index.")

    def transform(self, X, y=None):
        """Discretizes the data, ignoring any categorical features.

        Parameters
        ----------
        X : numeric array-like, shape (n_samples, n_features)
            Data containing features to be discretized.

        Returns
        -------
        T : numeric array-like, shape (n_samples, n_features)
            Categorical features will not be transformed.
        """
        check_is_fitted(self, ["min_", "ptp_", "clip_min_", "clip_max_"])
        X = check_array(X, dtype=float, ensure_2d=False)

        if X.ndim == 1:
            raise ValueError("Reshape your data.")

        if X.shape[1] != self.n_features_:
            raise ValueError("Transformed array does not have correct number "
                             "of features. Expecting {}, received {}."
                             .format(self.n_features_, X.shape[1]))

        with np.errstate(divide='ignore', invalid='ignore'):
            X_t = (X - self.min_) * self.n_bins_ // self.ptp_
            X_t[~np.isfinite(X_t)] = 0

        np.clip(X_t, self.clip_min_, self.clip_max_, out=X_t)
        return X_t
