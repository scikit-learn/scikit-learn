"""Transformers for missing value imputation"""
# Authors: Nicolas Tresegnie <nicolas.tresegnie@gmail.com>
#          Sergey Feldman <sergeyfeldman@gmail.com>
# License: BSD 3 clause

import warnings
import numbers

import numpy as np
import numpy.ma as ma
from scipy import sparse
from scipy import stats

from .base import BaseEstimator, TransformerMixin
from .utils import check_array
from .utils.sparsefuncs import _get_median
from .utils.validation import check_is_fitted
from .utils.validation import FLOAT_DTYPES
from .metrics import pairwise_distances
from .metrics.pairwise import _NAN_METRICS

from .neighbors.base import _check_weights
from .neighbors.base import _get_weights
from .utils.fixes import _object_dtype_isnan
from .utils import is_scalar_nan


__all__ = [
    'MissingIndicator',
    'SimpleImputer',
    'KNNImputer',
]


def _check_inputs_dtype(X, missing_values):
    if (X.dtype.kind in ("f", "i", "u") and
            not isinstance(missing_values, numbers.Real)):
        raise ValueError("'X' and 'missing_values' types are expected to be"
                         " both numerical. Got X.dtype={} and "
                         " type(missing_values)={}."
                         .format(X.dtype, type(missing_values)))


def _get_mask(X, value_to_mask):
    """Compute the boolean mask X == missing_values."""
    if is_scalar_nan(value_to_mask):
        if X.dtype.kind == "f":
            return np.isnan(X)
        elif X.dtype.kind in ("i", "u"):
            # can't have NaNs in integer array.
            return np.zeros(X.shape, dtype=bool)
        else:
            # np.isnan does not work on object dtypes.
            return _object_dtype_isnan(X)
    else:
        # X == value_to_mask with object dytpes does not always perform
        # element-wise for old versions of numpy
        return np.equal(X, value_to_mask)


def _most_frequent(array, extra_value, n_repeat):
    """Compute the most frequent value in a 1d array extended with
       [extra_value] * n_repeat, where extra_value is assumed to be not part
       of the array."""
    # Compute the most frequent value in array only
    if array.size > 0:
        with warnings.catch_warnings():
            # stats.mode raises a warning when input array contains objects due
            # to incapacity to detect NaNs. Irrelevant here since input array
            # has already been NaN-masked.
            warnings.simplefilter("ignore", RuntimeWarning)
            mode = stats.mode(array)

        most_frequent_value = mode[0][0]
        most_frequent_count = mode[1][0]
    else:
        most_frequent_value = 0
        most_frequent_count = 0

    # Compare to array + [extra_value] * n_repeat
    if most_frequent_count == 0 and n_repeat == 0:
        return np.nan
    elif most_frequent_count < n_repeat:
        return extra_value
    elif most_frequent_count > n_repeat:
        return most_frequent_value
    elif most_frequent_count == n_repeat:
        # Ties the breaks. Copy the behaviour of scipy.stats.mode
        if most_frequent_value < extra_value:
            return most_frequent_value
        else:
            return extra_value


class SimpleImputer(BaseEstimator, TransformerMixin):
    """Imputation transformer for completing missing values.

    Read more in the :ref:`User Guide <impute>`.

    Parameters
    ----------
    missing_values : number, string, np.nan (default) or None
        The placeholder for the missing values. All occurrences of
        `missing_values` will be imputed.

    strategy : string, optional (default="mean")
        The imputation strategy.

        - If "mean", then replace missing values using the mean along
          each column. Can only be used with numeric data.
        - If "median", then replace missing values using the median along
          each column. Can only be used with numeric data.
        - If "most_frequent", then replace missing using the most frequent
          value along each column. Can be used with strings or numeric data.
        - If "constant", then replace missing values with fill_value. Can be
          used with strings or numeric data.

        .. versionadded:: 0.20
           strategy="constant" for fixed value imputation.

    fill_value : string or numerical value, optional (default=None)
        When strategy == "constant", fill_value is used to replace all
        occurrences of missing_values.
        If left to the default, fill_value will be 0 when imputing numerical
        data and "missing_value" for strings or object data types.

    verbose : integer, optional (default=0)
        Controls the verbosity of the imputer.

    copy : boolean, optional (default=True)
        If True, a copy of X will be created. If False, imputation will
        be done in-place whenever possible. Note that, in the following cases,
        a new copy will always be made, even if `copy=False`:

        - If X is not an array of floating values;
        - If X is encoded as a CSR matrix.

    Attributes
    ----------
    statistics_ : array of shape (n_features,)
        The imputation fill value for each feature.

    Examples
    --------
    >>> import numpy as np
    >>> from sklearn.impute import SimpleImputer
    >>> imp_mean = SimpleImputer(missing_values=np.nan, strategy='mean')
    >>> imp_mean.fit([[7, 2, 3], [4, np.nan, 6], [10, 5, 9]])
    ... # doctest: +NORMALIZE_WHITESPACE
    SimpleImputer(copy=True, fill_value=None, missing_values=nan,
           strategy='mean', verbose=0)
    >>> X = [[np.nan, 2, 3], [4, np.nan, 6], [10, np.nan, 9]]
    >>> print(imp_mean.transform(X))
    ... # doctest: +NORMALIZE_WHITESPACE
    [[ 7.   2.   3. ]
     [ 4.   3.5  6. ]
     [10.   3.5  9. ]]

    Notes
    -----
    Columns which only contained missing values at `fit` are discarded upon
    `transform` if strategy is not "constant".

    """
    def __init__(self, missing_values=np.nan, strategy="mean",
                 fill_value=None, verbose=0, copy=True):
        self.missing_values = missing_values
        self.strategy = strategy
        self.fill_value = fill_value
        self.verbose = verbose
        self.copy = copy

    def _validate_input(self, X):
        allowed_strategies = ["mean", "median", "most_frequent", "constant"]
        if self.strategy not in allowed_strategies:
            raise ValueError("Can only use these strategies: {0} "
                             " got strategy={1}".format(allowed_strategies,
                                                        self.strategy))

        if self.strategy in ("most_frequent", "constant"):
            dtype = None
        else:
            dtype = FLOAT_DTYPES

        if not is_scalar_nan(self.missing_values):
            force_all_finite = True
        else:
            force_all_finite = "allow-nan"

        try:
            X = check_array(X, accept_sparse='csc', dtype=dtype,
                            force_all_finite=force_all_finite, copy=self.copy)
        except ValueError as ve:
            if "could not convert" in str(ve):
                raise ValueError("Cannot use {0} strategy with non-numeric "
                                 "data. Received datatype :{1}."
                                 "".format(self.strategy, X.dtype.kind))
            else:
                raise ve

        _check_inputs_dtype(X, self.missing_values)
        if X.dtype.kind not in ("i", "u", "f", "O"):
            raise ValueError("SimpleImputer does not support data with dtype "
                             "{0}. Please provide either a numeric array (with"
                             " a floating point or integer dtype) or "
                             "categorical data represented either as an array "
                             "with integer dtype or an array of string values "
                             "with an object dtype.".format(X.dtype))

        return X

    def fit(self, X, y=None):
        """Fit the imputer on X.

        Parameters
        ----------
        X : {array-like, sparse matrix}, shape (n_samples, n_features)
            Input data, where ``n_samples`` is the number of samples and
            ``n_features`` is the number of features.

        Returns
        -------
        self : SimpleImputer
        """
        X = self._validate_input(X)

        # default fill_value is 0 for numerical input and "missing_value"
        # otherwise
        if self.fill_value is None:
            if X.dtype.kind in ("i", "u", "f"):
                fill_value = 0
            else:
                fill_value = "missing_value"
        else:
            fill_value = self.fill_value

        # fill_value should be numerical in case of numerical input
        if (self.strategy == "constant" and
                X.dtype.kind in ("i", "u", "f") and
                not isinstance(fill_value, numbers.Real)):
            raise ValueError("'fill_value'={0} is invalid. Expected a "
                             "numerical value when imputing numerical "
                             "data".format(fill_value))

        if sparse.issparse(X):
            # missing_values = 0 not allowed with sparse data as it would
            # force densification
            if self.missing_values == 0:
                raise ValueError("Imputation not possible when missing_values "
                                 "== 0 and input is sparse. Provide a dense "
                                 "array instead.")
            else:
                self.statistics_ = self._sparse_fit(X,
                                                    self.strategy,
                                                    self.missing_values,
                                                    fill_value)
        else:
            self.statistics_ = self._dense_fit(X,
                                               self.strategy,
                                               self.missing_values,
                                               fill_value)

        return self

    def _sparse_fit(self, X, strategy, missing_values, fill_value):
        """Fit the transformer on sparse data."""
        mask_data = _get_mask(X.data, missing_values)
        n_implicit_zeros = X.shape[0] - np.diff(X.indptr)

        statistics = np.empty(X.shape[1])

        if strategy == "constant":
            # for constant strategy, self.statistcs_ is used to store
            # fill_value in each column
            statistics.fill(fill_value)

        else:
            for i in range(X.shape[1]):
                column = X.data[X.indptr[i]:X.indptr[i + 1]]
                mask_column = mask_data[X.indptr[i]:X.indptr[i + 1]]
                column = column[~mask_column]

                # combine explicit and implicit zeros
                mask_zeros = _get_mask(column, 0)
                column = column[~mask_zeros]
                n_explicit_zeros = mask_zeros.sum()
                n_zeros = n_implicit_zeros[i] + n_explicit_zeros

                if strategy == "mean":
                    s = column.size + n_zeros
                    statistics[i] = np.nan if s == 0 else column.sum() / s

                elif strategy == "median":
                    statistics[i] = _get_median(column,
                                                n_zeros)

                elif strategy == "most_frequent":
                    statistics[i] = _most_frequent(column,
                                                   0,
                                                   n_zeros)
        return statistics

    def _dense_fit(self, X, strategy, missing_values, fill_value):
        """Fit the transformer on dense data."""
        mask = _get_mask(X, missing_values)
        masked_X = ma.masked_array(X, mask=mask)

        # Mean
        if strategy == "mean":
            mean_masked = np.ma.mean(masked_X, axis=0)
            # Avoid the warning "Warning: converting a masked element to nan."
            mean = np.ma.getdata(mean_masked)
            mean[np.ma.getmask(mean_masked)] = np.nan

            return mean

        # Median
        elif strategy == "median":
            median_masked = np.ma.median(masked_X, axis=0)
            # Avoid the warning "Warning: converting a masked element to nan."
            median = np.ma.getdata(median_masked)
            median[np.ma.getmaskarray(median_masked)] = np.nan

            return median

        # Most frequent
        elif strategy == "most_frequent":
            # scipy.stats.mstats.mode cannot be used because it will no work
            # properly if the first element is masked and if its frequency
            # is equal to the frequency of the most frequent valid element
            # See https://github.com/scipy/scipy/issues/2636

            # To be able access the elements by columns
            X = X.transpose()
            mask = mask.transpose()

            if X.dtype.kind == "O":
                most_frequent = np.empty(X.shape[0], dtype=object)
            else:
                most_frequent = np.empty(X.shape[0])

            for i, (row, row_mask) in enumerate(zip(X[:], mask[:])):
                row_mask = np.logical_not(row_mask).astype(np.bool)
                row = row[row_mask]
                most_frequent[i] = _most_frequent(row, np.nan, 0)

            return most_frequent

        # Constant
        elif strategy == "constant":
            # for constant strategy, self.statistcs_ is used to store
            # fill_value in each column
            return np.full(X.shape[1], fill_value, dtype=X.dtype)

    def transform(self, X):
        """Impute all missing values in X.

        Parameters
        ----------
        X : {array-like, sparse matrix}, shape (n_samples, n_features)
            The input data to complete.
        """
        check_is_fitted(self, 'statistics_')

        X = self._validate_input(X)

        statistics = self.statistics_

        if X.shape[1] != statistics.shape[0]:
            raise ValueError("X has %d features per sample, expected %d"
                             % (X.shape[1], self.statistics_.shape[0]))

        # Delete the invalid columns if strategy is not constant
        if self.strategy == "constant":
            valid_statistics = statistics
        else:
            # same as np.isnan but also works for object dtypes
            invalid_mask = _get_mask(statistics, np.nan)
            valid_mask = np.logical_not(invalid_mask)
            valid_statistics = statistics[valid_mask]
            valid_statistics_indexes = np.flatnonzero(valid_mask)

            if invalid_mask.any():
                missing = np.arange(X.shape[1])[invalid_mask]
                if self.verbose:
                    warnings.warn("Deleting features without "
                                  "observed values: %s" % missing)
                X = X[:, valid_statistics_indexes]

        # Do actual imputation
        if sparse.issparse(X):
            if self.missing_values == 0:
                raise ValueError("Imputation not possible when missing_values "
                                 "== 0 and input is sparse. Provide a dense "
                                 "array instead.")
            else:
                mask = _get_mask(X.data, self.missing_values)
                indexes = np.repeat(np.arange(len(X.indptr) - 1, dtype=np.int),
                                    np.diff(X.indptr))[mask]

                X.data[mask] = valid_statistics[indexes].astype(X.dtype,
                                                                copy=False)
        else:
            mask = _get_mask(X, self.missing_values)
            n_missing = np.sum(mask, axis=0)
            values = np.repeat(valid_statistics, n_missing)
            coordinates = np.where(mask.transpose())[::-1]

            X[coordinates] = values

        return X


class MissingIndicator(BaseEstimator, TransformerMixin):
    """Binary indicators for missing values.

    Note that this component typically should not not be used in a vanilla
    :class:`Pipeline` consisting of transformers and a classifier, but rather
    could be added using a :class:`FeatureUnion` or :class:`ColumnTransformer`.

    Read more in the :ref:`User Guide <impute>`.

    Parameters
    ----------
    missing_values : number, string, np.nan (default) or None
        The placeholder for the missing values. All occurrences of
        `missing_values` will be indicated (True in the output array), the
        other values will be marked as False.

    features : str, optional
        Whether the imputer mask should represent all or a subset of
        features.

        - If "missing-only" (default), the imputer mask will only represent
          features containing missing values during fit time.
        - If "all", the imputer mask will represent all features.

    sparse : boolean or "auto", optional
        Whether the imputer mask format should be sparse or dense.

        - If "auto" (default), the imputer mask will be of same type as
          input.
        - If True, the imputer mask will be a sparse matrix.
        - If False, the imputer mask will be a numpy array.

    error_on_new : boolean, optional
        If True (default), transform will raise an error when there are
        features with missing values in transform that have no missing values
        in fit. This is applicable only when ``features="missing-only"``.

    Attributes
    ----------
    features_ : ndarray, shape (n_missing_features,) or (n_features,)
        The features indices which will be returned when calling ``transform``.
        They are computed during ``fit``. For ``features='all'``, it is
        to ``range(n_features)``.

    Examples
    --------
    >>> import numpy as np
    >>> from sklearn.impute import MissingIndicator
    >>> X1 = np.array([[np.nan, 1, 3],
    ...                [4, 0, np.nan],
    ...                [8, 1, 0]])
    >>> X2 = np.array([[5, 1, np.nan],
    ...                [np.nan, 2, 3],
    ...                [2, 4, 0]])
    >>> indicator = MissingIndicator()
    >>> indicator.fit(X1)  # doctest: +NORMALIZE_WHITESPACE
    MissingIndicator(error_on_new=True, features='missing-only',
             missing_values=nan, sparse='auto')
    >>> X2_tr = indicator.transform(X2)
    >>> X2_tr
    array([[False,  True],
           [ True, False],
           [False, False]])

    """

    def __init__(self, missing_values=np.nan, features="missing-only",
                 sparse="auto", error_on_new=True):
        self.missing_values = missing_values
        self.features = features
        self.sparse = sparse
        self.error_on_new = error_on_new

    def _get_missing_features_info(self, X):
        """Compute the imputer mask and the indices of the features
        containing missing values.

        Parameters
        ----------
        X : {ndarray or sparse matrix}, shape (n_samples, n_features)
            The input data with missing values. Note that ``X`` has been
            checked in ``fit`` and ``transform`` before to call this function.

        Returns
        -------
        imputer_mask : {ndarray or sparse matrix}, shape \
(n_samples, n_features) or (n_samples, n_features_with_missing)
            The imputer mask of the original data.

        features_with_missing : ndarray, shape (n_features_with_missing)
            The features containing missing values.

        """
        if sparse.issparse(X) and self.missing_values != 0:
            mask = _get_mask(X.data, self.missing_values)

            # The imputer mask will be constructed with the same sparse format
            # as X.
            sparse_constructor = (sparse.csr_matrix if X.format == 'csr'
                                  else sparse.csc_matrix)
            imputer_mask = sparse_constructor(
                (mask, X.indices.copy(), X.indptr.copy()),
                shape=X.shape, dtype=bool)

            missing_values_mask = imputer_mask.copy()
            missing_values_mask.eliminate_zeros()
            features_with_missing = (
                np.flatnonzero(np.diff(missing_values_mask.indptr))
                if missing_values_mask.format == 'csc'
                else np.unique(missing_values_mask.indices))

            if self.sparse is False:
                imputer_mask = imputer_mask.toarray()
            elif imputer_mask.format == 'csr':
                imputer_mask = imputer_mask.tocsc()
        else:
            if sparse.issparse(X):
                # case of sparse matrix with 0 as missing values. Implicit and
                # explicit zeros are considered as missing values.
                X = X.toarray()
            imputer_mask = _get_mask(X, self.missing_values)
            features_with_missing = np.flatnonzero(imputer_mask.sum(axis=0))

            if self.sparse is True:
                imputer_mask = sparse.csc_matrix(imputer_mask)

        return imputer_mask, features_with_missing

    def fit(self, X, y=None):
        """Fit the transformer on X.

        Parameters
        ----------
        X : {array-like, sparse matrix}, shape (n_samples, n_features)
            Input data, where ``n_samples`` is the number of samples and
            ``n_features`` is the number of features.

        Returns
        -------
        self : object
            Returns self.
        """
        if not is_scalar_nan(self.missing_values):
            force_all_finite = True
        else:
            force_all_finite = "allow-nan"
        X = check_array(X, accept_sparse=('csc', 'csr'),
                        force_all_finite=force_all_finite)
        _check_inputs_dtype(X, self.missing_values)

        self._n_features = X.shape[1]

        if self.features not in ('missing-only', 'all'):
            raise ValueError("'features' has to be either 'missing-only' or "
                             "'all'. Got {} instead.".format(self.features))

        if not ((isinstance(self.sparse, str) and
                self.sparse == "auto") or isinstance(self.sparse, bool)):
            raise ValueError("'sparse' has to be a boolean or 'auto'. "
                             "Got {!r} instead.".format(self.sparse))

        self.features_ = (self._get_missing_features_info(X)[1]
                          if self.features == 'missing-only'
                          else np.arange(self._n_features))

        return self

    def transform(self, X):
        """Generate missing values indicator for X.

        Parameters
        ----------
        X : {array-like, sparse matrix}, shape (n_samples, n_features)
            The input data to complete.

        Returns
        -------
        Xt : {ndarray or sparse matrix}, shape (n_samples, n_features)
            The missing indicator for input data. The data type of ``Xt``
            will be boolean.

        """
        check_is_fitted(self, "features_")

        if not is_scalar_nan(self.missing_values):
            force_all_finite = True
        else:
            force_all_finite = "allow-nan"
        X = check_array(X, accept_sparse=('csc', 'csr'),
                        force_all_finite=force_all_finite)
        _check_inputs_dtype(X, self.missing_values)

        if X.shape[1] != self._n_features:
            raise ValueError("X has a different number of features "
                             "than during fitting.")

        imputer_mask, features = self._get_missing_features_info(X)

        if self.features == "missing-only":
            features_diff_fit_trans = np.setdiff1d(features, self.features_)
            if (self.error_on_new and features_diff_fit_trans.size > 0):
                raise ValueError("The features {} have missing values "
                                 "in transform but have no missing values "
                                 "in fit.".format(features_diff_fit_trans))

            if (self.features_.size > 0 and
                    self.features_.size < self._n_features):
                imputer_mask = imputer_mask[:, self.features_]

        return imputer_mask

    def fit_transform(self, X, y=None):
        """Generate missing values indicator for X.

        Parameters
        ----------
        X : {array-like, sparse matrix}, shape (n_samples, n_features)
            The input data to complete.

        Returns
        -------
        Xt : {ndarray or sparse matrix}, shape (n_samples, n_features)
            The missing indicator for input data. The data type of ``Xt``
            will be boolean.

        """
        return self.fit(X, y).transform(X)


class KNNImputer(BaseEstimator, TransformerMixin):
    """Imputation for completing missing values using k-Nearest Neighbors.

    Each sample's missing values are imputed using values from ``n_neighbors``
    nearest neighbors found in the training set. Each missing feature is then
    imputed as the average, either weighted or unweighted, of these neighbors.
    Note that if a sample has more than one feature missing, then the
    neighbors for that sample can be different depending on the particular
    feature being imputed. Finally, where the number of donor neighbors is
    less than ``n_neighbors``, the training set average for that feature is
    used during imputation.

    Parameters
    ----------
    missing_values : number, string, np.nan (default) or None
        The placeholder for the missing values. All occurrences of
        `missing_values` will be imputed.

    n_neighbors : int, optional (default = 5)
        Number of neighboring samples to use for imputation.

    weights : str or callable, optional (default = "uniform")
        Weight function used in prediction.  Possible values:

        - 'uniform' : uniform weights.  All points in each neighborhood
          are weighted equally.
        - 'distance' : weight points by the inverse of their distance.
          in this case, closer neighbors of a query point will have a
          greater influence than neighbors which are further away.
        - [callable] : a user-defined function which accepts an
          array of distances, and returns an array of the same shape
          containing the weights.

    metric : str or callable, optional (default = "masked_euclidean")
        Distance metric for searching neighbors. Possible values:
        - 'masked_euclidean'
        - [callable] : a user-defined function which conforms to the
          definition of _pairwise_callable(X, Y, metric, **kwds). In other
          words, the function accepts two arrays, X and Y, and a
          ``missing_values`` keyword in **kwds and returns a scalar distance
          value.

    row_max_missing : float, optional (default = 0.5)
        The maximum fraction of columns (i.e. features) that can be missing
        before the sample is excluded from nearest neighbor imputation. It
        means that such rows will not be considered a potential donor in
        ``fit()``, and in ``transform()`` their missing feature values will be
        imputed to be the column mean for the entire dataset.

    col_max_missing : float, optional (default = 0.8)
        The maximum fraction of rows (or samples) that can be missing
        for any feature beyond which an error is raised.

    copy : boolean, optional (default = True)
        If True, a copy of X will be created. If False, imputation will
        be done in-place whenever possible. Note that, if metric is
        "masked_euclidean" and copy=False then missing_values in the
        input matrix X will be overwritten with zeros.

    Attributes
    ----------
    statistics_ : 1-D array of length {n_features}
        The 1-D array contains the mean of each feature calculated using
        observed (i.e. non-missing) values. This is used for imputing
        missing values in samples that are either excluded from nearest
        neighbors search because they have too many ( > row_max_missing)
        missing features or because all of the sample's k-nearest neighbors
        (i.e., the potential donors) also have the relevant feature value
        missing.

    References
    ----------
    * Olga Troyanskaya, Michael Cantor, Gavin Sherlock, Pat Brown, Trevor
      Hastie, Robert Tibshirani, David Botstein and Russ B. Altman, Missing
      value estimation methods for DNA microarrays, BIOINFORMATICS Vol. 17
      no. 6, 2001 Pages 520-525.

    Examples
    --------
    >>> from sklearn.impute import KNNImputer
    >>> nan = float("NaN")
    >>> X = [[1, 2, nan], [3, 4, 3], [nan, 6, 5], [8, 8, 7]]
    >>> imputer = KNNImputer(n_neighbors=2, weights="uniform")
    >>> imputer.fit_transform(X)
    array([[1. , 2. , 4. ],
           [3. , 4. , 3. ],
           [5.5, 6. , 5. ],
           [8. , 8. , 7. ]])
    """

    def __init__(self, missing_values=np.nan, n_neighbors=5,
                 weights="uniform", metric="masked_euclidean",
                 row_max_missing=0.5, col_max_missing=0.8, copy=True):

        self.missing_values = missing_values
        self.n_neighbors = n_neighbors
        self.weights = weights
        self.metric = metric
        self.row_max_missing = row_max_missing
        self.col_max_missing = col_max_missing
        self.copy = copy

    def _impute(self, dist, X, fitted_X, mask, mask_fx):
        """Helper function to find and impute missing values"""

        # For each column, find and impute
        n_rows_X, n_cols_X = X.shape
        for c in range(n_cols_X):
            if not np.any(mask[:, c], axis=0):
                continue

            # Row index for receivers and potential donors (pdonors)
            receivers_row_idx = np.where(mask[:, c])[0]
            pdonors_row_idx = np.where(~mask_fx[:, c])[0]

            # Impute using column mean if n_neighbors are not available
            if len(pdonors_row_idx) < self.n_neighbors:
                warnings.warn("Insufficient number of neighbors! "
                              "Filling in column mean.")
                X[receivers_row_idx, c] = self.statistics_[c]
                continue

            # Get distance from potential donors
            dist_pdonors = dist[receivers_row_idx][:, pdonors_row_idx]
            dist_pdonors = dist_pdonors.reshape(-1,
                                                len(pdonors_row_idx))

            # Argpartition to separate actual donors from the rest
            pdonors_idx = np.argpartition(
                dist_pdonors, self.n_neighbors - 1, axis=1)

            # Get final donors row index from pdonors
            donors_idx = pdonors_idx[:, :self.n_neighbors]

            # Get weights or None
            dist_pdonors_rows = np.arange(len(donors_idx))[:, None]
            weight_matrix = _get_weights(
                dist_pdonors[
                    dist_pdonors_rows, donors_idx], self.weights)
            donor_row_idx_ravel = donors_idx.ravel()

            # Retrieve donor values and calculate kNN score
            fitted_X_temp = fitted_X[pdonors_row_idx]
            donors = fitted_X_temp[donor_row_idx_ravel, c].reshape(
                (-1, self.n_neighbors))
            donors_mask = _get_mask(donors, self.missing_values)
            donors = np.ma.array(donors, mask=donors_mask)

            # Final imputation
            imputed = np.ma.average(donors, axis=1,
                                    weights=weight_matrix)
            X[receivers_row_idx, c] = imputed.data
        return X

    def fit(self, X, y=None):
        """Fit the imputer on X.

        Parameters
        ----------
        X : array-like, shape (n_samples, n_features)
            Input data, where ``n_samples`` is the number of samples and
            ``n_features`` is the number of features.

        Returns
        -------
        self : object
            Returns self.
        """

        # Check data integrity and calling arguments
        if not is_scalar_nan(self.missing_values):
            force_all_finite = True
        else:
            force_all_finite = "allow-nan"
            if self.metric not in _NAN_METRICS and not callable(
                    self.metric):
                raise ValueError(
                    "The selected metric does not support NaN values.")
        X = check_array(X, accept_sparse=False, dtype=np.float64,
                        force_all_finite=force_all_finite, copy=self.copy)
        self.weights = _check_weights(self.weights)

        # Check if % missing in any column > col_max_missing
        mask = _get_mask(X, self.missing_values)
        if np.any(mask.sum(axis=0) > (X.shape[0] * self.col_max_missing)):
            raise ValueError("Some column(s) have more than {}% missing values"
                             .format(self.col_max_missing * 100))
        X_col_means = np.ma.array(X, mask=mask).mean(axis=0).data

        # Check if % missing in any row > row_max_missing
        bad_rows = mask.sum(axis=1) > (mask.shape[1] * self.row_max_missing)
        if np.any(bad_rows):
            warnings.warn(
                "There are rows with more than {0}% missing values. These "
                "rows are not included as donor neighbors."
                    .format(self.row_max_missing * 100))

            # Remove rows that have more than row_max_missing % missing
            X = X[~bad_rows, :]

        # Check if sufficient neighboring samples available
        if X.shape[0] < self.n_neighbors:
            raise ValueError("There are only %d samples, but n_neighbors=%d."
                             % (X.shape[0], self.n_neighbors))
        self.fitted_X_ = X
        self.statistics_ = X_col_means

        return self

    def transform(self, X):
        """Impute all missing values in X.

        Parameters
        ----------
        X : array-like, shape = [n_samples, n_features]
            The input data to complete.

        Returns
        -------
        X : array-like, shape = [n_samples, n_features]
            The imputed dataset.
        """

        check_is_fitted(self, ["fitted_X_", "statistics_"])
        if not is_scalar_nan(self.missing_values):
            force_all_finite = True
        else:
            force_all_finite = "allow-nan"
        X = check_array(X, accept_sparse=False, dtype=FLOAT_DTYPES,
                        force_all_finite=force_all_finite, copy=self.copy)

        # Get fitted data and ensure correct dimension
        n_rows_fit_X, n_cols_fit_X = self.fitted_X_.shape
        n_rows_X, n_cols_X = X.shape

        if n_cols_X != n_cols_fit_X:
            raise ValueError("Incompatible dimension between the fitted "
                             "dataset and the one to be transformed.")
        mask = _get_mask(X, self.missing_values)

        row_total_missing = mask.sum(axis=1)
        if not np.any(row_total_missing):
            return X

        # Check for excessive missingness in rows
        bad_rows = row_total_missing > (mask.shape[1] * self.row_max_missing)
        if np.any(bad_rows):
            warnings.warn(
                "There are rows with more than {0}% missing values. The "
                "missing features in these rows are imputed with column means."
                    .format(self.row_max_missing * 100))
            X_bad = X[bad_rows, :]
            X = X[~bad_rows, :]
            mask = mask[~bad_rows]
            row_total_missing = mask.sum(axis=1)
        row_has_missing = row_total_missing.astype(np.bool)

        if np.any(row_has_missing):

            # Mask for fitted_X
            mask_fx = _get_mask(self.fitted_X_, self.missing_values)

            # Pairwise distances between receivers and fitted samples
            dist = np.empty((len(X), len(self.fitted_X_)))
            dist[row_has_missing] = pairwise_distances(
                X[row_has_missing], self.fitted_X_, metric=self.metric,
                squared=False, missing_values=self.missing_values)

            # Find and impute missing
            X = self._impute(dist, X, self.fitted_X_, mask, mask_fx)

        # Merge bad rows to X and mean impute their missing values
        if np.any(bad_rows):
            bad_missing_index = np.where(_get_mask(X_bad, self.missing_values))
            X_bad[bad_missing_index] = np.take(self.statistics_,
                                               bad_missing_index[1])
            X_merged = np.empty((n_rows_X, n_cols_X))
            X_merged[bad_rows, :] = X_bad
            X_merged[~bad_rows, :] = X
            X = X_merged
        return X

    def fit_transform(self, X, y=None, **fit_params):
        """Fit KNNImputer and impute all missing values in X.

        Parameters
        ----------
        X : array-like, shape (n_samples, n_features)
            Input data, where ``n_samples`` is the number of samples and
            ``n_features`` is the number of features.

        Returns
        -------
        X : array-like, shape (n_samples, n_features)
            Returns imputed dataset.
        """
        return self.fit(X).transform(X)
