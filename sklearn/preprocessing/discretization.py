# Author: Henry Lin <hlin117@gmail.com>

from __future__ import division, absolute_import

import numbers
import numpy as np
import warnings

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

    ignored_features : int array-like (default=None)
        Column indices of ignored features. (Example: Categorical features.)
        If None, all features will be discretized.

    Attributes
    ----------
    min_ : float array, shape (n_features_,)
        Minimum value per feature in X. An ignored feature at index i will
        have min_[i] == 0.

    ptp_ : float array, shape (n_features_,)
        X.max(axis=0) - X.min(axis=0). An ignored feature at index i will
        have ptp_[i] == 1.

    n_bins_ : int array, shape (n_features_,)
        Number of bins per feature. An ignored feature at index i will
        have n_bins_[i] == 1.

    Example
    -------
    >>> X = [[-2, 1, -4,   -1],
    ...      [-1, 2, -3, -0.5],
    ...      [ 0, 3, -2,  0.5],
    ...      [ 1, 4, -1,    2]]
    >>> est = KBinsDiscretizer(n_bins=3)
    >>> est.fit(X) # doctest: +ELLIPSIS, +NORMALIZE_WHITESPACE
    KBinsDiscretizer(...)
    >>> est.transform(X) # doctest: +SKIP
    array([[ 0., 0., 0., 0.],
           [ 1., 1., 1., 0.],
           [ 2., 2., 2., 1.],
           [ 2., 2., 2., 2.]])
    """

    def __init__(self, n_bins=2, ignored_features=None):
        self.n_bins = n_bins
        self.ignored_features = ignored_features

    def fit(self, X, y=None):
        """Fits the estimator.

        Parameters
        ----------
        X : numeric array-like, shape (n_samples, n_features)
            Data to be discretized.

        y : ignored

        Returns
        -------
        self
        """
        X = check_array(X, dtype=float)

        n_features = X.shape[1]
        ignored_features = self._check_ignored_features(self.ignored_features,
                                                        X.shape[1])

        # Manipulating min and max so categorical data is not affected by
        # transformation.
        min = np.min(X, axis=0)
        max = np.max(X, axis=0)

        min[ignored_features] = 0
        max[ignored_features] = 1

        self.min_ = min
        self.ptp_ = max - min

        n_bins = self._check_n_bins(self.n_bins, n_features)
        n_bins[ignored_features] = 1
        self.n_bins_ = n_bins

        # Clipping features ensures outliers are within boundary after
        # transformation, but also ensures categorical features can be ignored.
        numeric_features = np.delete(np.arange(n_features), ignored_features)
        clip_min = np.repeat(-np.inf, n_features)
        clip_min[numeric_features] = 0

        clip_max = np.repeat(np.inf, n_features)
        clip_max[numeric_features] = self.n_bins_[numeric_features] - 1

        self._clip_min = clip_min
        self._clip_max = clip_max

        same_min_max = np.where(self.ptp_ == 0)[0]
        if len(same_min_max) > 0:
            warnings.warn("Fitted X contained constant features at indices "
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
            return np.ones(n_features) * n_bins

        n_bins = check_array(n_bins, dtype=int, copy=True, ensure_2d=False)

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
        return n_bins

    @staticmethod
    def _check_ignored_features(ignored_features, n_features):
        if ignored_features is None:
            return np.array([], dtype='int64')

        ignored_features = check_array(ignored_features, ensure_2d=False,
                                       dtype=int)
        ignored_features = column_or_1d(ignored_features)

        if len(set(ignored_features)) != ignored_features.shape[0]:
            raise ValueError("Duplicate ignored column indices found.")

        if np.all(ignored_features >= 0) and \
                np.all(ignored_features < n_features):
            return ignored_features

        raise ValueError("Invalid ignored feature index.")

    def transform(self, X, y=None):
        """Discretizes the data.

        Parameters
        ----------
        X : numeric array-like, shape (n_samples, n_features)
            Data containing features to be discretized.

        Returns
        -------
        Xt : numeric array-like, shape (n_samples, n_features)
        """
        check_is_fitted(self, ["min_", "ptp_", "_clip_min", "_clip_max"])
        X = check_array(X, dtype=float)

        n_features = self.n_bins_.shape[0]
        if X.shape[1] != n_features:
            raise ValueError("Incorrect number of features. Expecting {}, "
                             "received {}.".format(n_features, X.shape[1]))

        with np.errstate(divide='ignore', invalid='ignore'):
            Xt = (X - self.min_) * self.n_bins_ // self.ptp_
            Xt[~np.isfinite(Xt)] = 0

        np.clip(Xt, self._clip_min, self._clip_max, out=Xt)
        return Xt
