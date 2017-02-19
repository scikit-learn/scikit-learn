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
        Column indexes of categorical features. Categorical features
        will not be discretized and will be ignored by the transformation.

    Attributes
    ----------
    mins_ : float array, shape (n_features_,)
        Minimum value per feature in X.

    maxes_ : float array, shape (n_features_,)
        Maximum value per feature in X.

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
    >>> dis.transform(X) # doctest: +SKIP
    array([[0, 0, 0, 0],
           [1, 1, 1, 0],
           [2, 2, 2, 1],
           [2, 2, 2, 2]])

    Categorical features can be specified in the constructor.

    >>> dis = KBinsDiscretizer(n_bins=3, categorical_features=[1])
    >>> dis.fit(X) # doctest: +ELLIPSIS, +NORMALIZE_WHITESPACE
    KBinsDiscretizer(...)

    Discretization will ignore categorical features.
    >>> dis.transform(X) # doctest: +SKIP
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

        self._check_categorical_features(self.categorical_features, X.shape[1])

        n_features = X.shape[1]

        if isinstance(self.n_bins, numbers.Number) and self.n_bins < 2:
            raise ValueError("Discretizer received an invalid number "
                             "of bins. Received {0}, expected at least 2."
                             .format(self.n_bins))
        else:
            check_array(self.n_bins, dtype=int)

        self.mins_ = np.min(X, axis=0)
        self.maxes_ = np.max(X, axis=0)
        self.n_features_ = n_features
        return self

    @staticmethod
    def _check_categorical_features(cat_features, n_features):
        if cat_features is None:
            return

        cat_features = check_array(cat_features, ensure_2d=False, dtype=int)
        cat_features = column_or_1d(cat_features)

        if len(set(cat_features)) != cat_features.shape[0]:
            raise ValueError("Duplicate categorical column indexes found.")

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
        T : float array-like, shape (n_samples, n_features)
            Categorical features will not be transformed.
        """
        check_is_fitted(self, ["mins_", "maxes_", "n_features_"])
        X = check_array(X, dtype=float, ensure_2d=False)

        if X.ndim == 1:
            raise ValueError("Reshape your data.")

        if X.shape[1] != self.n_features_:
            raise ValueError("Transformed array does not have correct number "
                             "of features. Expecting {}, received {}."
                             .format(self.n_features_, X.shape[1]))

        if self.categorical_features:
            cat_features = np.asarray(sorted(self.categorical_features))
        else:
            cat_features = np.array([], dtype='int64')

        # Manipulate the input attributes so categorical features can be
        # ignored in transformation process. Also avoid divide by zero errors.
        same_min_max = np.arange(self.n_features_)[self.maxes_ == self.mins_]
        if len(same_min_max) != 0:
            import warnings
            warnings.warn("Fitted X contained continuous features at indexes "
                          "{}. These features will be discretized to the "
                          "same value."
                          .format(", ".join(str(i) for i in same_min_max)))

        ignored_features = np.concatenate((cat_features, same_min_max))

        mins = self.mins_.copy()
        mins[ignored_features] = 0

        maxes = self.maxes_.copy()
        maxes[ignored_features] = 1

        if isinstance(self.n_bins, numbers.Number):
            n_bins = np.ones(self.n_features_) * self.n_bins
        else:
            n_bins = self.n_bins.copy()
        n_bins[ignored_features] = 1

        X_t = (X - mins) * n_bins // (maxes - mins)

        # Need to make sure outliers are within boundary, but also make sure
        # categorical features can be ignored.
        numeric_features = np.delete(np.arange(self.n_features_), cat_features)
        clip_min = np.repeat(-np.inf, self.n_features_)
        clip_min[numeric_features] = 0

        clip_max = np.repeat(np.inf, self.n_features_)
        clip_max[numeric_features] = self.n_bins - 1

        np.clip(X_t, clip_min, clip_max, out=X_t)

        X_t[:, same_min_max] = 0
        return X_t
