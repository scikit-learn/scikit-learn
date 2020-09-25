# Authors: Henry Lin <hlin117@gmail.com>
#          Tom Dupré la Tour
#          Juan Carlos Alfaro Jiménez <JuanCarlos.Alfaro@uclm.es>

# License: BSD


import numbers
import numpy as np
import warnings

from . import OneHotEncoder

from ..base import BaseEstimator, TransformerMixin
from ._mdlp_discretization import find_binning_thresholds
from ._mdlp_discretization import DTYPE, ITYPE, SIZE
from ..utils.multiclass import check_classification_targets
from ..utils.validation import check_array
from ..utils.validation import check_is_fitted
from ..utils.validation import _deprecate_positional_args


class _BaseDiscretizer(TransformerMixin, BaseEstimator):
    """Base class for all discretizers."""

    def _validate_dtype(self, X):
        """Validate the output data type."""
        supported_dtypes = (None, np.float32, np.float64)

        if self.dtype is None:
            output_dtype = X.dtype
        elif self.dtype in supported_dtypes:
            output_dtype = self.dtype
        else:
            raise ValueError("Valid options for 'dtype' are "
                             f"type are {supported_dtypes}. "
                             f"Got dtype={self.dtype}  instead.")

        return output_dtype

    def _validate_encode(self):
        """Validate the method to encode the result."""
        valid_encodes = ("onehot", "onehot-dense", "ordinal")

        if self.encode not in valid_encodes:
            raise ValueError("Valid options for 'encode' "
                             f"are {valid_encodes}. Got "
                             f"encode='{self.encode}' instead.")

    def transform(self, X):
        """Discretize the data.

        Parameters
        ----------
        X : array-like of shape (n_samples, n_features), dtype={int, float}
            The data to discretize.

        Returns
        -------
        Xt : {ndarray, sparse matrix}, dtype={np.float32, np.float64}
            The data in the binned space. Will be a sparse matrix
            if `self.encode="onehot"` and ndarray otherwise.
        """
        check_is_fitted(self)

        # Check the input data and the attributes
        dtype = (np.float64, np.float32) if self.dtype is None else self.dtype

        Xt = self._validate_data(X, reset=False, copy=True, dtype=dtype)

        bin_edges = self.bin_edges_

        for jj in range(self.n_features_in_):
            # Values close to a bin edge are susceptible to numeric
            # instability. Thus, add epsilon to ensure these values
            # are binned correctly w.r.t. to the decimal truncation
            rtol = 1.e-5
            atol = 1.e-8
            eps = atol + rtol * np.abs(Xt[:, jj])
            Xt[:, jj] = np.digitize(Xt[:, jj] + eps, bin_edges[jj][1:])

        np.clip(Xt, 0, self.n_bins_ - 1, out=Xt)

        if self.encode == "ordinal":
            return Xt

        dtype_init = None

        if "onehot" in self.encode:
            dtype_init = self._encoder.dtype
            self._encoder.dtype = Xt.dtype
        try:
            Xt_enc = self._encoder.transform(Xt)
        finally:
            # Revert the initial type to avoid modifying
            self._encoder.dtype = dtype_init

        return Xt_enc

    def inverse_transform(self, Xt):
        """Transform discretized data back to original feature space.

        Note that this function does not regenerate the
        original data due to discretization rounding.

        Parameters
        ----------
        Xt : array-like of shape (n_samples, n_features), dtype={int, float}
            The transformed data in the binned space.

        Returns
        -------
        Xinv : ndarray, dtype={np.float32, np.float64}
            The data in the original feature space.
        """
        check_is_fitted(self)

        if "onehot" in self.encode:
            Xt = self._encoder.inverse_transform(Xt)

        Xinv = self._validate_data(Xt,
                                   copy=True,
                                   dtype=(np.float64, np.float32))

        for jj in range(self.n_features_in_):
            bin_edges = self.bin_edges_[jj]
            bin_centers = (bin_edges[1:] + bin_edges[:-1]) * 0.5
            Xinv[:, jj] = bin_centers[np.int_(Xinv[:, jj])]

        return Xinv


class KBinsDiscretizer(_BaseDiscretizer):
    """Bin continuous data into intervals using a given number of bins.

    Read more in the :ref:`User Guide <preprocessing_discretization>`.

    .. versionadded:: 0.20

    Parameters
    ----------
    n_bins : int or array-like of shape (n_features,), default=5
        The number of bins to produce. Raises ValueError if ``n_bins < 2``.

    encode : {'onehot', 'onehot-dense', 'ordinal'}, default='onehot'
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

    strategy : {'uniform', 'quantile', 'kmeans'}, default='quantile'
        Strategy used to define the widths of the bins.

        uniform
            All bins in each feature have identical widths.
        quantile
            All bins in each feature have the same number of points.
        kmeans
            Values in each bin have the same nearest center of a 1D k-means
            cluster.

    dtype : {np.float32, np.float64}, default=None
        The desired data-type for the output. If None, output dtype is
        consistent with input dtype. Only np.float32 and np.float64 are
        supported.

        .. versionadded:: 0.24

    Attributes
    ----------
    n_bins_ : ndarray of shape (n_features,), dtype=np.int_
        Number of bins per feature. Bins whose width are too small
        (i.e., <= 1e-8) are removed with a warning.

    bin_edges_ : ndarray of ndarray of shape (n_features,)
        The edges of each bin. Contain arrays of varying shapes ``(n_bins_, )``
        Ignored features will have empty arrays.

    See Also
    --------
    Binarizer : Class to bin values as `0` or `1` based on a threshold.
    MDLPDiscretizer : Class to bin values using the MDLP strategy.

    Notes
    -----
    In bin edges for feature ``i``, the first and last values are used only for
    ``inverse_transform``. During transform, bin edges are extended to::

      np.concatenate([-np.inf, bin_edges_[i][1:-1], np.inf])

    You can combine ``KBinsDiscretizer`` with
    :class:`~sklearn.compose.ColumnTransformer` if you only want to preprocess
    part of the features.

    ``KBinsDiscretizer`` might produce constant features (e.g., when
    ``encode = 'onehot'`` and certain bins do not contain any data).
    These features can be removed with feature selection algorithms
    (e.g., :class:`~sklearn.feature_selection.VarianceThreshold`).

    Examples
    --------
    >>> X = [[-2, 1, -4,   -1],
    ...      [-1, 2, -3, -0.5],
    ...      [ 0, 3, -2,  0.5],
    ...      [ 1, 4, -1,    2]]
    >>> est = KBinsDiscretizer(n_bins=3, encode='ordinal', strategy='uniform')
    >>> est.fit(X)
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
    """

    @_deprecate_positional_args
    def __init__(self, n_bins=5, *, encode='onehot', strategy='quantile',
                 dtype=None):
        self.n_bins = n_bins
        self.encode = encode
        self.strategy = strategy
        self.dtype = dtype

    def fit(self, X, y=None):
        """
        Fit the estimator.

        Parameters
        ----------
        X : array-like of shape (n_samples, n_features), dtype={int, float}
            Data to be discretized.

        y : None
            Ignored. This parameter exists only for compatibility with
            :class:`~sklearn.pipeline.Pipeline`.

        Returns
        -------
        self
        """
        X = self._validate_data(X, dtype='numeric')

        output_dtype = self._validate_dtype(X)
        self._validate_encode()

        valid_strategy = ('uniform', 'quantile', 'kmeans')
        if self.strategy not in valid_strategy:
            raise ValueError("Valid options for 'strategy' are {}. "
                             "Got strategy={!r} instead."
                             .format(valid_strategy, self.strategy))

        n_features = X.shape[1]
        n_bins = self._validate_n_bins(n_features)

        bin_edges = np.zeros(n_features, dtype=object)
        for jj in range(n_features):
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
                bin_edges[jj] = np.asarray(np.percentile(column, quantiles))

            elif self.strategy == 'kmeans':
                from ..cluster import KMeans  # fixes import loops

                # Deterministic initialization with uniform spacing
                uniform_edges = np.linspace(col_min, col_max, n_bins[jj] + 1)
                init = (uniform_edges[1:] + uniform_edges[:-1])[:, None] * 0.5

                # 1D k-means procedure
                km = KMeans(n_clusters=n_bins[jj], init=init, n_init=1)
                centers = km.fit(column[:, None]).cluster_centers_[:, 0]
                # Must sort, centers may be unsorted even with sorted init
                centers.sort()
                bin_edges[jj] = (centers[1:] + centers[:-1]) * 0.5
                bin_edges[jj] = np.r_[col_min, bin_edges[jj], col_max]

            # Remove bins whose width are too small (i.e., <= 1e-8)
            if self.strategy in ('quantile', 'kmeans'):
                mask = np.ediff1d(bin_edges[jj], to_begin=np.inf) > 1e-8
                bin_edges[jj] = bin_edges[jj][mask]
                if len(bin_edges[jj]) - 1 != n_bins[jj]:
                    warnings.warn('Bins whose width are too small (i.e., <= '
                                  '1e-8) in feature %d are removed. Consider '
                                  'decreasing the number of bins.' % jj)
                    n_bins[jj] = len(bin_edges[jj]) - 1

        self.bin_edges_ = bin_edges
        self.n_bins_ = n_bins

        if 'onehot' in self.encode:
            self._encoder = OneHotEncoder(
                categories=[np.arange(i) for i in self.n_bins_],
                sparse=self.encode == 'onehot',
                dtype=output_dtype)
            # Fit the OneHotEncoder with toy datasets
            # so that it's ready for use after the KBinsDiscretizer is fitted
            self._encoder.fit(np.zeros((1, len(self.n_bins_))))

        return self

    def _validate_n_bins(self, n_features):
        """Returns n_bins_, the number of bins per feature.
        """
        orig_bins = self.n_bins
        if isinstance(orig_bins, numbers.Number):
            if not isinstance(orig_bins, numbers.Integral):
                raise ValueError("{} received an invalid n_bins type. "
                                 "Received {}, expected int."
                                 .format(KBinsDiscretizer.__name__,
                                         type(orig_bins).__name__))
            if orig_bins < 2:
                raise ValueError("{} received an invalid number "
                                 "of bins. Received {}, expected at least 2."
                                 .format(KBinsDiscretizer.__name__, orig_bins))
            return np.full(n_features, orig_bins, dtype=int)

        n_bins = check_array(orig_bins, dtype=int, copy=True,
                             ensure_2d=False)

        if n_bins.ndim > 1 or n_bins.shape[0] != n_features:
            raise ValueError("n_bins must be a scalar or array "
                             "of shape (n_features,).")

        bad_nbins_value = (n_bins < 2) | (n_bins != orig_bins)

        violating_indices = np.where(bad_nbins_value)[0]
        if violating_indices.shape[0] > 0:
            indices = ", ".join(str(i) for i in violating_indices)
            raise ValueError("{} received an invalid number "
                             "of bins at indices {}. Number of bins "
                             "must be at least 2, and must be an int."
                             .format(KBinsDiscretizer.__name__, indices))
        return n_bins


class MDLPDiscretizer(_BaseDiscretizer):
    """Bin continuous data into intervals using the MDLP strategy.

    Read more in the :ref:`User Guide <preprocessing_discretization>`.

    .. note::

      This estimator is still **experimental** for now: the predictions
      and the API might change without any deprecation cycle. To use it,
      you need to explicitly import ``enable_mdlp_discretizer``::

        >>> # Explictly require this experimental feature
        >>> from sklearn.experimental import enable_mdlp_discretizer  # noqa
        >>> # Now you can import normally from preprocessing
        >>> from sklearn.preprocessing import MDLPDiscretizer

    Parameters
    ----------
    encode : {"onehot", "onehot-dense", "ordinal"}, default="onehot"
        The method used to encode the transformed result.

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

    dtype : {np.float32, np.float64}, default=None
        The desired data type for the output. If `None`, the output data
        type is consistent with input data type. Only `np.float32` and
        `np.float64` are supported.

    Attributes
    ----------
    n_bins_ : ndarray of shape (n_features,), dtype=int
        The number of bins per feature.

    bin_edges_ : list of ndarray of shape (n_features,) of dtype=float
        The edges of each bin, containing arrays of varying
        shapes, according to the attribute `n_bins_`.

    See Also
    --------
    Binarizer : Class to bin values as `0` or `1` based on a threshold.
    KBinsDiscretizer : Class to bin values using a given number of intervals.

    Notes
    -----
    In bin edges for feature `i`, the first and last values are used only
    for `inverse_transform`. During `transform`, they are extended to::

        np.concatenate([-np.inf, bin_edges_[i][1:-1], np.inf])

    If you only want to preprocess part of the features, you can combine
    `MDLPDiscretizer` with :class:`~sklearn.compose.ColumnTransformer`.

    Examples
    --------
    >>> X = [[-2, 1, -4,   -1],
    ...      [-1, 2, -3, -0.5],
    ...      [ 0, 3, -2,  0.5],
    ...      [ 1, 4, -1,    2]]
    >>> y = [0, 0, 0, 1]
    >>> discretizer = MDLPDiscretizer(encode="ordinal")
    >>> discretizer.fit(X, y)
    MDLPDiscretizer(...)
    >>> Xt = discretizer.transform(X)
    >>> Xt  # doctest: +SKIP
    array([[ 0., 0., 0., 0.],
           [ 0., 0., 0., 0.],
           [ 0., 0., 0., 0.],
           [ 1., 1., 1., 1.]])

    Sometimes is useful to convert the data back into the original feature
    space. The `inverse_transform` function converts the binned data into
    the original feature space. Each value will be equal to the mean of
    the two bin edges.

    >>> discretizer.bin_edges_[0]
    array([-2. ,  0.5,  1. ])
    >>> discretizer.inverse_transform(Xt)
    array([[-0.75 ,  2.25 , -2.75 ,  0.125],
           [-0.75 ,  2.25 , -2.75 ,  0.125],
           [-0.75 ,  2.25 , -2.75 ,  0.125],
           [ 0.75 ,  3.75 , -1.25 ,  1.625]])
    """

    def __init__(self, *, encode="onehot", dtype=None):
        self.encode = encode
        self.dtype = dtype

    def fit(self, X, y):
        """Fit the discretizer to the training data and labels.

        Parameters
        ----------
        X : array-like of shape (n_samples, n_features)
            The data to determine the bins for each feature.

        y : array-like of shape (n_samples,) or (n_samples, n_outputs)
            The labels to determine the bins for each feature.

        Returns
        -------
        self : MDLPDiscretizer
            The fitted discretizer to the training data and labels.
        """
        X, y = self._validate_data(X, y, multi_output=True)

        # The discretizer is available only for classification
        check_classification_targets(y)

        output_dtype = self._validate_dtype(X)
        self._validate_encode()

        if y.ndim == 1:
            # Reshape is necessary to preserve the
            # data contiguity against "np.newaxis"
            y = np.reshape(y, (-1, 1))

        # Determine output settings
        n_samples, n_outputs = y.shape
        n_features = self.n_features_in_

        n_classes = np.zeros(n_outputs, dtype=SIZE)
        y_encoded = np.zeros((n_samples, n_outputs), dtype=ITYPE)

        for k in range(n_outputs):
            classes_k, y_encoded[:, k] = np.unique(
                y[:, k], return_inverse=True)

            n_classes[k] = classes_k.shape[0]

        X = np.ascontiguousarray(X, dtype=DTYPE)
        X_idx_sorted = np.argsort(X, axis=0)
        y_encoded = np.ascontiguousarray(y_encoded, dtype=ITYPE)

        bin_edges, n_bins = find_binning_thresholds(X,
                                                    X_idx_sorted,
                                                    y_encoded,
                                                    n_outputs,
                                                    n_classes)

        self.n_bins_ = n_bins
        self.bin_edges_ = np.zeros(n_features, dtype=object)

        for j in range(n_features):
            col_min, col_max = np.min(X[:, j]), np.max(X[:, j])

            if col_min == col_max or self.n_bins_[j] == 1:
                if col_min == col_max:
                    warnings.warn(f"Feature {j} is constant and "
                                  "will be replaced with 0.")
                else:
                    warnings.warn("No bin edges found for feature "
                                  f"{j}. It will be replaced with 0.")

                self.n_bins_[j] = 1
                self.bin_edges_[j] = np.array([-np.inf, np.inf])

                continue

            self.bin_edges_[j] = np.concatenate([[col_min],
                                                 bin_edges[j],
                                                 [col_max]])

        # Fit an encoder with a "toy" dataset
        # to use after the estimator is fitted
        if "onehot" in self.encode:
            categories = [np.arange(i) for i in self.n_bins_]
            sparse = self.encode == "onehot"
            X = np.zeros((1, self.n_bins_.shape[0]))

            self._encoder = OneHotEncoder(sparse=sparse,
                                          categories=categories,
                                          dtype=output_dtype)

            self._encoder.fit(X)

        return self
