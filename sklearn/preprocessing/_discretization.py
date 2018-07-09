# -*- coding: utf-8 -*-

# Author: Henry Lin <hlin117@gmail.com>
#         Tom Dupré la Tour

# License: BSD

from __future__ import division, absolute_import

import numbers
import numpy as np
import warnings

from . import OneHotEncoder
from ._encoders import _transform_selected

from ..base import BaseEstimator, TransformerMixin
from ..utils.validation import check_array
from ..utils.validation import check_is_fitted
from ..utils.validation import column_or_1d
from ..utils.fixes import np_version


class KBinsDiscretizer(BaseEstimator, TransformerMixin):
    """Bin continuous data into intervals.

    Read more in the :ref:`User Guide <preprocessing_discretization>`.

    Parameters
    ----------
    n_bins : int or array-like, shape (n_features,) (default=5)
        The number of bins to produce. The intervals for the bins are
        determined by the minimum and maximum of the input data.
        Raises ValueError if ``n_bins < 2``.

        If ``n_bins`` is an array, and there is an ignored feature at
        index ``i``, ``n_bins[i]`` will be ignored.

    ignored_features : int array-like (default=None)
        Column indices of ignored features. (Example: Categorical features.)
        If ``None``, all features will be discretized.

    encode : {'onehot', 'onehot-dense', 'ordinal'}, (default='onehot')
        Method used to encode the transformed result.

        onehot
            Encode the transformed result with one-hot encoding
            and return a sparse matrix. Ignored features are always
            stacked to the right.
        onehot-dense
            Encode the transformed result with one-hot encoding
            and return a dense array. Ignored features are always
            stacked to the right.
        ordinal
            Return the bin identifier encoded as an integer value.

    dtype : number type, default=np.float
        Desired dtype of output.

    strategy : {'uniform', 'quantile', 'kmeans'}, (default='quantile')
        Strategy used to define the widths of the bins.

        uniform
            All bins in each feature have identical widths.
        quantile
            All bins in each feature have the same number of points.
        kmeans
            Values in each bin have the same nearest center of a 1D k-means
            cluster.

    Attributes
    ----------
    n_bins_ : int array, shape (n_features,)
        Number of bins per feature. An ignored feature at index ``i``
        will have ``n_bins_[i] == 0``.

    bin_edges_ : array of arrays, shape (n_features, )
        The edges of each bin. Contain arrays of varying shapes (n_bins_, ).
        Ignored features will have empty arrays.

    transformed_features_ : int array, shape (n_features,)
        Features which are transformed.

    Examples
    --------
    >>> X = [[-2, 1, -4,   -1],
    ...      [-1, 2, -3, -0.5],
    ...      [ 0, 3, -2,  0.5],
    ...      [ 1, 4, -1,    2]]
    >>> est = KBinsDiscretizer(n_bins=3, encode='ordinal', strategy='uniform')
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
    data into the original feature space. Each value will be equal to the mean
    of the two bin edges.

    >>> est.bin_edges_[0]
    array([-2., -1.,  0.,  1.])
    >>> est.inverse_transform(Xt)
    array([[-1.5,  1.5, -3.5, -0.5],
           [-0.5,  2.5, -2.5, -0.5],
           [ 0.5,  3.5, -1.5,  0.5],
           [ 0.5,  3.5, -1.5,  1.5]])

    Notes
    -----
    In bin edges for feature ``i``, the first and last values are used only for
    ``inverse_transform``. During transform, bin edges are extended to::

      np.concatenate([-np.inf, bin_edges_[i][1:-1], np.inf])

    See also
    --------
     sklearn.preprocessing.Binarizer : class used to bin values as ``0`` or
        ``1`` based on a parameter ``threshold``.
    """

    def __init__(self, n_bins=5, ignored_features=None, encode='onehot',
                 dtype=np.float64, strategy='quantile'):
        self.n_bins = n_bins
        self.ignored_features = ignored_features
        self.encode = encode
        self.dtype = dtype
        self.strategy = strategy

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

        valid_encode = ('onehot', 'onehot-dense', 'ordinal')
        if self.encode not in valid_encode:
            raise ValueError("Valid options for 'encode' are {}. "
                             "Got encode={!r} instead."
                             .format(valid_encode, self.encode))
        valid_strategy = ('uniform', 'quantile', 'kmeans')
        if self.strategy not in valid_strategy:
            raise ValueError("Valid options for 'strategy' are {}. "
                             "Got strategy={!r} instead."
                             .format(valid_strategy, self.strategy))

        n_features = X.shape[1]
        ignored = self._validate_ignored_features(n_features)
        self.transformed_features_ = np.delete(np.arange(n_features), ignored)

        n_bins = self._validate_n_bins(n_features, ignored)

        bin_edges = np.zeros(n_features, dtype=object)
        for jj in range(n_features):
            if jj in ignored:
                bin_edges[jj] = np.array([])
                continue
            column = X[:, jj]
            col_min, col_max = column.min(), column.max()

            if col_min == col_max:
                warnings.warn("Feature %d is constant and will be "
                              "replaced with 0." % jj)
                n_bins[jj] = 1
                bin_edges[jj] = np.array([-np.inf, np.inf])
                continue

            if self.strategy == 'uniform':
                bin_edges[jj] = np.linspace(col_min, col_max, n_bins[jj] + 1)

            elif self.strategy == 'quantile':
                quantiles = np.linspace(0, 100, n_bins[jj] + 1)
                if np_version < (1, 9):
                    quantiles = list(quantiles)
                bin_edges[jj] = np.asarray(np.percentile(column, quantiles))

            elif self.strategy == 'kmeans':
                from ..cluster import KMeans  # fixes import loops

                # Deterministic initialization with uniform spacing
                uniform_edges = np.linspace(col_min, col_max, n_bins[jj] + 1)
                init = (uniform_edges[1:] + uniform_edges[:-1])[:, None] * 0.5

                # 1D k-means procedure
                km = KMeans(n_clusters=n_bins[jj], init=init, n_init=1)
                centers = km.fit(column[:, None]).cluster_centers_[:, 0]
                bin_edges[jj] = (centers[1:] + centers[:-1]) * 0.5
                bin_edges[jj] = np.r_[col_min, bin_edges[jj], col_max]

        self.bin_edges_ = bin_edges
        self.n_bins_ = n_bins

        return self

    def _validate_n_bins(self, n_features, ignored):
        """Returns n_bins_, the number of bins per feature.

        Also ensures that ignored bins are zero.
        """
        orig_bins = self.n_bins
        if isinstance(orig_bins, numbers.Number):
            if not isinstance(orig_bins, (numbers.Integral, np.integer)):
                raise ValueError("{} received an invalid n_bins type. "
                                 "Received {}, expected int."
                                 .format(KBinsDiscretizer.__name__,
                                         type(orig_bins).__name__))
            if orig_bins < 2:
                raise ValueError("{} received an invalid number "
                                 "of bins. Received {}, expected at least 2."
                                 .format(KBinsDiscretizer.__name__, orig_bins))
            return np.ones(n_features, dtype=np.int) * orig_bins

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
        Xt : numeric array-like or sparse matrix
            Data in the binned space.
        """
        check_is_fitted(self, ["bin_edges_"])
        X = self._validate_X_post_fit(X)

        Xt = _transform_selected(X, self._transform, self.dtype,
                                 self.transformed_features_, copy=True,
                                 retain_order=True)

        if self.encode == 'ordinal':
            return Xt

        # Only one-hot encode discretized features
        mask = np.ones(X.shape[1], dtype=bool)
        if self.ignored_features is not None:
            mask[self.ignored_features] = False

        encode_sparse = self.encode == 'onehot'
        return OneHotEncoder(n_values=self.n_bins_[mask],
                             categorical_features='all'
                             if self.ignored_features is None else mask,
                             sparse=encode_sparse).fit_transform(Xt)

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

        bin_edges = self.bin_edges_[trans]
        for jj in range(X.shape[1]):
            # Values which are close to a bin edge are susceptible to numeric
            # instability. Add eps to X so these values are binned correctly
            # with respect to their decimal truncation. See documentation of
            # numpy.isclose for an explanation of ``rtol`` and ``atol``.
            rtol = 1.e-5
            atol = 1.e-8
            eps = atol + rtol * np.abs(X[:, jj])

            X[:, jj] = np.digitize(X[:, jj] + eps, bin_edges[jj][1:])

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
        check_is_fitted(self, ["bin_edges_"])

        # Currently, OneHotEncoder doesn't support inverse_transform
        if self.encode != 'ordinal':
            raise ValueError("inverse_transform only supports "
                             "'encode = ordinal'. Got encode={!r} instead."
                             .format(self.encode))

        Xt = self._validate_X_post_fit(Xt)
        trans = self.transformed_features_
        Xinv = Xt.copy()
        Xinv_sel = Xinv[:, trans]

        n_features = Xinv_sel.shape[1]
        for jj in range(n_features):
            bin_edges = self.bin_edges_[trans][jj]
            bin_centers = (bin_edges[1:] + bin_edges[:-1]) * 0.5
            Xinv_sel[:, jj] = bin_centers[np.int_(Xinv_sel[:, jj])]

        Xinv[:, trans] = Xinv_sel
        return Xinv
