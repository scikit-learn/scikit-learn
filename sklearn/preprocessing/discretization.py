# Author: Henry Lin <hlin117@gmail.com>

# License: BSD

from __future__ import division, absolute_import

import numbers
import numpy as np
import warnings

from sklearn.base import BaseEstimator, TransformerMixin
from sklearn.preprocessing.data import _transform_selected
from sklearn.utils.validation import check_array
from sklearn.utils.validation import check_is_fitted
from sklearn.utils.validation import column_or_1d


class KBinsDiscretizer(BaseEstimator, TransformerMixin):
    """Bins continuous data into k equal width intervals.

    Parameters
    ----------
    n_bins : int or array-like, shape (n_features,) (default=2)
        The number of bins to produce. The intervals for the bins are
        determined by the minimum and maximum of the input data.
        Raises ValueError if ``n_bins < 2``.

        If ``n_bins`` is an array, and there is an ignored feature at
        index ``i``, ``n_bins[i]`` will be ignored.

    ignored_features : int array-like (default=None)
        Column indices of ignored features. (Example: Categorical features.)
        If ``None``, all features will be discretized.

    Attributes
    ----------
    offset_ : float array, shape (n_features,)
        Minimum value per feature in ``X``. An ignored feature at index
        ``i`` will have ``offset_[i] == 0``.

    n_bins_ : int array, shape (n_features,)
        Number of bins per feature. An ignored feature at index ``i``
        will have ``n_bins_[i] == 0``.

    bin_width_ : float array, shape (n_features,)
        The width of each bin. Ignored features will have widths equal
        to ``0``.

    transformed_features_ : int array, shape (n_features,)
        Features which are transformed.

    Example
    -------
    >>> X = [[-2, 1, -4,   -1],
    ...      [-1, 2, -3, -0.5],
    ...      [ 0, 3, -2,  0.5],
    ...      [ 1, 4, -1,    2]]
    >>> est = KBinsDiscretizer(n_bins=3)
    >>> est.fit(X)  # doctest: +ELLIPSIS, +NORMALIZE_WHITESPACE
    KBinsDiscretizer(...)
    >>> Xt = est.transform(X)
    >>> Xt  # doctest: +SKIP
    array([[ 0., 0., 0., 0.],
           [ 1., 1., 1., 0.],
           [ 2., 2., 2., 1.],
           [ 2., 2., 2., 2.]])

    Sometimes it may be useful to convert the data back into the original
    feature space. The ``inverse_transform`` function converts the binned
    data into the original feature space. Each value will be a distance
    of ``bin_width_ / 2`` from the nearest bin edge.

    >>> est.bin_width_
    array([ 1.,  1.,  1.,  1.])
    >>> est.inverse_transform(Xt)
    array([[-1.5,  1.5, -3.5, -0.5],
           [-0.5,  2.5, -2.5, -0.5],
           [ 0.5,  3.5, -1.5,  0.5],
           [ 0.5,  3.5, -1.5,  1.5]])

    Notes
    -----
    Bin edges for feature ``i`` are defined as::

      np.concatenate([
        -np.inf,
        offset_[i] + bin_width_[i] * np.arange(1, n_bins_[i]),
        np.inf
      ])

    See also
    --------
     sklearn.preprocessing.Binarizer : class used to bin values as ``0`` or
        ``1`` based on a parameter ``threshold``.
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
        X = check_array(X, dtype='numeric')

        n_features = X.shape[1]
        ignored = self._validate_ignored_features(n_features)
        self.transformed_features_ = np.delete(np.arange(n_features), ignored)

        offset = np.min(X, axis=0)
        offset[ignored] = 0
        self.offset_ = offset
        self.n_bins_ = self._validate_n_bins(n_features, ignored)

        ptp = np.ptp(X, axis=0)
        same_min_max = np.where(ptp == 0)[0]
        if len(same_min_max) > 0:
            warnings.warn("Features {} are constant and will be replaced with "
                          "0.".format(", ".join(str(i) for i in same_min_max)))

        with np.errstate(divide='ignore', invalid='ignore'):
            bin_widths = ptp / self.n_bins_
        bin_widths[ignored] = 0
        self.bin_width_ = bin_widths
        return self

    def _validate_n_bins(self, n_features, ignored):
        """Returns n_bins_, the number of bins per feature.

        Also ensures that ignored bins are zero.
        """
        orig_bins = self.n_bins
        if isinstance(orig_bins, numbers.Number):
            if not isinstance(orig_bins, np.int):
                raise ValueError("{} received an invalid n_bins type. "
                                 "Received {}, expected int."
                                 .format(KBinsDiscretizer.__name__,
                                         type(orig_bins).__name__))
            if orig_bins < 2:
                raise ValueError("{} received an invalid number "
                                 "of bins. Received {}, expected at least 2."
                                 .format(KBinsDiscretizer.__name__, orig_bins))
            return np.ones(n_features) * orig_bins

        n_bins = check_array(orig_bins, dtype=np.int, copy=True,
                             ensure_2d=False)

        if n_bins.ndim > 1 or n_bins.shape[0] != n_features:
            raise ValueError("n_bins must be a scalar or array "
                             "of shape (n_features,).")

        bad_nbins_value = (n_bins < 2) | (n_bins != orig_bins)
        bad_nbins_value[ignored] = False

        violating_indices = np.where(bad_nbins_value)[0]
        if violating_indices.shape[0] > 0:
            indices = ", ".join(str(i) for i in violating_indices)
            raise ValueError("{} received an invalid number "
                             "of bins at indices {}. Number of bins "
                             "must be at least 2, and must be an int."
                             .format(KBinsDiscretizer.__name__, indices))
        n_bins[ignored] = 0
        return n_bins

    def _validate_ignored_features(self, n_features):
        ignored = self.ignored_features
        if ignored is None:
            return np.array([], dtype='int64')

        ignored = check_array(ignored, ensure_2d=False, dtype=int)
        ignored = column_or_1d(ignored)

        if len(set(ignored)) != ignored.shape[0]:
            raise ValueError("Duplicate ignored column indices found.")

        if np.all(ignored >= 0) and np.all(ignored < n_features):
            return ignored

        raise ValueError("Invalid ignored feature index.")

    def transform(self, X):
        """Discretizes the data.

        Parameters
        ----------
        X : numeric array-like, shape (n_samples, n_features)
            Data to be discretized.

        Returns
        -------
        Xt : numeric array-like, shape (n_samples, n_features)
            Data in the binned space.
        """
        check_is_fitted(self, ["offset_", "bin_width_"])
        X = self._validate_X_post_fit(X)

        return _transform_selected(X, self._transform,
                                   self.transformed_features_, copy=True,
                                   retain_order=True)

    def _validate_X_post_fit(self, X):
        X = check_array(X, dtype='numeric')

        n_features = self.n_bins_.shape[0]
        if X.shape[1] != n_features:
            raise ValueError("Incorrect number of features. Expecting {}, "
                             "received {}.".format(n_features, X.shape[1]))
        return X

    def _transform(self, X):
        """Performs transformation on X, with no ignored features."""
        trans = self.transformed_features_

        X -= self.offset_[trans]
        bin_width = self.bin_width_[trans]

        # Rescale into [-1, bin_width] range
        with np.errstate(divide='ignore', invalid='ignore'):
            X /= bin_width

        # Values which are a multiple of the bin width are susceptible to
        # numeric instability. Add eps to X so these values are binned
        # correctly. See documentation for numpy.isclose for an explanation
        # of ``rtol`` and ``atol``.
        rtol = 1.e-5
        atol = 1.e-8
        eps = atol + rtol * bin_width
        np.floor(X + eps, out=X)

        X[~np.isfinite(X)] = 0  # Case when a feature is constant
        np.clip(X, 0, self.n_bins_[trans] - 1, out=X)
        return X

    def inverse_transform(self, Xt):
        """Transforms discretized data back to original feature space.

        Note that this function does not regenerate the original data
        due to discretization rounding.

        Parameters
        ----------
        Xt : numeric array-like, shape (n_sample, n_features)
            Transformed data in the binned space.

        Returns
        -------
        Xinv : numeric array-like
            Data in the original feature space.
        """
        check_is_fitted(self, ["offset_", "bin_width_"])
        Xt = self._validate_X_post_fit(Xt)
        trans = self.transformed_features_
        Xinv = Xt.copy()
        Xinv_sel = Xinv[:, trans]

        Xinv_sel += 0.5
        Xinv_sel *= self.bin_width_[trans]
        Xinv_sel += self.offset_[trans]

        Xinv[:, trans] = Xinv_sel
        return Xinv
