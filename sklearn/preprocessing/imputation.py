# Authors: Nicolas Tresegnie <nicolas.tresegnie@gmail.com>
#          Ashim Bhattarai <"ashimb9" + "\100gmail\56com">
# License: BSD 3 clause
from __future__ import division
import warnings

import numpy as np
import numpy.ma as ma
from scipy import sparse
from scipy import stats

from ..base import BaseEstimator, TransformerMixin
from ..utils import check_array
from ..utils.sparsefuncs import _get_median
from ..utils.validation import check_is_fitted
from ..utils.validation import FLOAT_DTYPES

from ..externals import six

zip = six.moves.zip
map = six.moves.map

__all__ = [
    'Imputer',
    'KNNImputer'
]


def _get_mask(X, value_to_mask):
    """Compute the boolean mask X == missing_values."""
    if value_to_mask == "NaN" or np.isnan(value_to_mask):
        return np.isnan(X)
    else:
        return X == value_to_mask


def _most_frequent(array, extra_value, n_repeat):
    """Compute the most frequent value in a 1d array extended with
       [extra_value] * n_repeat, where extra_value is assumed to be not part
       of the array."""
    # Compute the most frequent value in array only
    if array.size > 0:
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


# Code for function _unique1d taken directly from Numpy
def _unique1d(ar, return_index=False, return_inverse=False,
              return_counts=False):
    """
    Find the unique elements of an array, ignoring shape.
    """
    ar = np.asanyarray(ar).flatten()

    optional_indices = return_index or return_inverse
    optional_returns = optional_indices or return_counts

    if ar.size == 0:
        if not optional_returns:
            ret = ar
        else:
            ret = (ar,)
            if return_index:
                ret += (np.empty(0, np.bool),)
            if return_inverse:
                ret += (np.empty(0, np.bool),)
            if return_counts:
                ret += (np.empty(0, np.intp),)
        return ret

    if optional_indices:
        perm = ar.argsort(kind='mergesort' if return_index else 'quicksort')
        aux = ar[perm]
    else:
        ar.sort()
        aux = ar
    flag = np.concatenate(([True], aux[1:] != aux[:-1]))

    if not optional_returns:
        ret = aux[flag]
    else:
        ret = (aux[flag],)
        if return_index:
            ret += (perm[flag],)
        if return_inverse:
            iflag = np.cumsum(flag) - 1
            inv_idx = np.empty(ar.shape, dtype=np.intp)
            inv_idx[perm] = iflag
            ret += (inv_idx,)
        if return_counts:
            idx = np.concatenate(np.nonzero(flag) + ([ar.size],))
            ret += (np.diff(idx),)
    return ret


# Code for function _unique taken directly from Numpy
def _unique(ar, return_index=False, return_inverse=False,
            return_counts=False, axis=None):
    """
    Find the unique elements of an array.

    Returns the sorted unique elements of an array. There are three optional
    outputs in addition to the unique elements: the indices of the input array
    that give the unique values, the indices of the unique array that
    reconstruct the input array, and the number of times each unique value
    comes up in the input array.

    Parameters
    ----------
    ar : array_like
        Input array. Unless `axis` is specified, this will be flattened if it
        is not already 1-D.
    return_index : bool, optional
        If True, also return the indices of `ar` (along the specified axis,
        if provided, or in the flattened array) that result in the unique
        array.
    return_inverse : bool, optional
        If True, also return the indices of the unique array (for the specified
        axis, if provided) that can be used to reconstruct `ar`.
    return_counts : bool, optional
        If True, also return the number of times each unique item appears
        in `ar`.
        .. versionadded:: 1.9.0
    axis : int or None, optional
        The axis to operate on. If None, `ar` will be flattened beforehand.
        Otherwise, duplicate items will be removed along the provided axis,
        with all the other axes belonging to the each of the unique elements.
        Object arrays or structured arrays that contain objects are not
        supported if the `axis` kwarg is used.
        .. versionadded:: 1.13.0



    Returns
    -------
    unique : ndarray
        The sorted unique values.
    unique_indices : ndarray, optional
        The indices of the first occurrences of the unique values in the
        original array. Only provided if `return_index` is True.
    unique_inverse : ndarray, optional
        The indices to reconstruct the original array from the
        unique array. Only provided if `return_inverse` is True.
    unique_counts : ndarray, optional
        The number of times each of the unique values comes up in the
        original array. Only provided if `return_counts` is True.
        .. versionadded:: 1.9.0

    See Also
    --------
    numpy.lib.arraysetops : Module with a number of other functions for
                            performing set operations on arrays.

    Examples
    --------
    >>> np.unique([1, 1, 2, 2, 3, 3])
    array([1, 2, 3])
    >>> a = np.array([[1, 1], [2, 3]])
    >>> np.unique(a)
    array([1, 2, 3])

    Return the unique rows of a 2D array

    >>> a = np.array([[1, 0, 0], [1, 0, 0], [2, 3, 4]])
    >>> np.unique(a, axis=0)
    array([[1, 0, 0], [2, 3, 4]])

    Return the indices of the original array that give the unique values:

    >>> a = np.array(['a', 'b', 'b', 'c', 'a'])
    >>> u, indices = np.unique(a, return_index=True)
    >>> u
    array(['a', 'b', 'c'],
           dtype='|S1')
    >>> indices
    array([0, 1, 3])
    >>> a[indices]
    array(['a', 'b', 'c'],
           dtype='|S1')

    Reconstruct the input array from the unique values:

    >>> a = np.array([1, 2, 6, 4, 2, 3, 2])
    >>> u, indices = np.unique(a, return_inverse=True)
    >>> u
    array([1, 2, 3, 4, 6])
    >>> indices
    array([0, 1, 4, 3, 1, 2, 1])
    >>> u[indices]
    array([1, 2, 6, 4, 2, 3, 2])

    """
    ar = np.asanyarray(ar)
    if axis is None:
        return _unique1d(ar, return_index, return_inverse, return_counts)
    if not (-ar.ndim <= axis < ar.ndim):
        raise ValueError('Invalid axis kwarg specified for unique')

    ar = np.swapaxes(ar, axis, 0)
    orig_shape, orig_dtype = ar.shape, ar.dtype
    # Must reshape to a contiguous 2D array for this to work...
    ar = ar.reshape(orig_shape[0], -1)
    ar = np.ascontiguousarray(ar)

    if ar.dtype.char in (np.typecodes['AllInteger'] +
                         np.typecodes['Datetime'] + 'S'):
        # Optimization: Creating a view of your data with a np.void data type
        # of size the number of bytes in a full row. Handles any type where
        # items have a unique binary representation, i.e. 0 is only 0,
        # not +0 and -0.
        dtype = np.dtype((np.void, ar.dtype.itemsize * ar.shape[1]))
    else:
        dtype = [('f{i}'.format(i=i), ar.dtype) for i in range(ar.shape[1])]

    try:
        consolidated = ar.view(dtype)
    except TypeError:
        # There's no good way to do this for object arrays, etc...
        msg = 'The axis argument to unique is not supported for dtype {dt}'
        raise TypeError(msg.format(dt=ar.dtype))

    def reshape_uniq(uniq):
        uniq = uniq.view(orig_dtype)
        uniq = uniq.reshape(-1, *orig_shape[1:])
        uniq = np.swapaxes(uniq, 0, axis)
        return uniq

    output = _unique1d(consolidated, return_index,
                       return_inverse, return_counts)
    if not (return_index or return_inverse or return_counts):
        return reshape_uniq(output)
    else:
        uniq = reshape_uniq(output[0])
        return (uniq,) + output[1:]


class Imputer(BaseEstimator, TransformerMixin):
    """Imputation transformer for completing missing values.

    Read more in the :ref:`User Guide <imputation>`.

    Parameters
    ----------
    missing_values : integer or "NaN", optional (default="NaN")
        The placeholder for the missing values. All occurrences of
        `missing_values` will be imputed. For missing values encoded as np.nan,
        use the string value "NaN".

    strategy : string, optional (default="mean")
        The imputation strategy.

        - If "mean", then replace missing values using the mean along
          the axis.
        - If "median", then replace missing values using the median along
          the axis.
        - If "most_frequent", then replace missing using the most frequent
          value along the axis.

    axis : integer, optional (default=0)
        The axis along which to impute.

        - If `axis=0`, then impute along columns.
        - If `axis=1`, then impute along rows.

    verbose : integer, optional (default=0)
        Controls the verbosity of the imputer.

    copy : boolean, optional (default=True)
        If True, a copy of X will be created. If False, imputation will
        be done in-place whenever possible. Note that, in the following cases,
        a new copy will always be made, even if `copy=False`:

        - If X is not an array of floating values;
        - If X is sparse and `missing_values=0`;
        - If `axis=0` and X is encoded as a CSR matrix;
        - If `axis=1` and X is encoded as a CSC matrix.

    Attributes
    ----------
    statistics_ : array of shape (n_features,)
        The imputation fill value for each feature if axis == 0.

    Notes
    -----
    - When ``axis=0``, columns which only contained missing values at `fit`
      are discarded upon `transform`.
    - When ``axis=1``, an exception is raised if there are rows for which it is
      not possible to fill in the missing values (e.g., because they only
      contain missing values).
    """
    def __init__(self, missing_values="NaN", strategy="mean",
                 axis=0, verbose=0, copy=True):
        self.missing_values = missing_values
        self.strategy = strategy
        self.axis = axis
        self.verbose = verbose
        self.copy = copy

    def fit(self, X, y=None):
        """Fit the imputer on X.

        Parameters
        ----------
        X : {array-like, sparse matrix}, shape (n_samples, n_features)
            Input data, where ``n_samples`` is the number of samples and
            ``n_features`` is the number of features.

        Returns
        -------
        self : Imputer
            Returns self.
        """
        # Check parameters
        allowed_strategies = ["mean", "median", "most_frequent"]
        if self.strategy not in allowed_strategies:
            raise ValueError("Can only use these strategies: {0} "
                             " got strategy={1}".format(allowed_strategies,
                                                        self.strategy))

        if self.axis not in [0, 1]:
            raise ValueError("Can only impute missing values on axis 0 and 1, "
                             " got axis={0}".format(self.axis))

        # Since two different arrays can be provided in fit(X) and
        # transform(X), the imputation data will be computed in transform()
        # when the imputation is done per sample (i.e., when axis=1).
        if self.axis == 0:
            X = check_array(X, accept_sparse='csc', dtype=np.float64,
                            force_all_finite=False)

            if sparse.issparse(X):
                self.statistics_ = self._sparse_fit(X,
                                                    self.strategy,
                                                    self.missing_values,
                                                    self.axis)
            else:
                self.statistics_ = self._dense_fit(X,
                                                   self.strategy,
                                                   self.missing_values,
                                                   self.axis)

        return self

    def _sparse_fit(self, X, strategy, missing_values, axis):
        """Fit the transformer on sparse data."""
        # Imputation is done "by column", so if we want to do it
        # by row we only need to convert the matrix to csr format.
        if axis == 1:
            X = X.tocsr()
        else:
            X = X.tocsc()

        # Count the zeros
        if missing_values == 0:
            n_zeros_axis = np.zeros(X.shape[not axis], dtype=int)
        else:
            n_zeros_axis = X.shape[axis] - np.diff(X.indptr)

        # Mean
        if strategy == "mean":
            if missing_values != 0:
                n_non_missing = n_zeros_axis

                # Mask the missing elements
                mask_missing_values = _get_mask(X.data, missing_values)
                mask_valids = np.logical_not(mask_missing_values)

                # Sum only the valid elements
                new_data = X.data.copy()
                new_data[mask_missing_values] = 0
                X = sparse.csc_matrix((new_data, X.indices, X.indptr),
                                      copy=False)
                sums = X.sum(axis=0)

                # Count the elements != 0
                mask_non_zeros = sparse.csc_matrix(
                    (mask_valids.astype(np.float64),
                     X.indices,
                     X.indptr), copy=False)
                s = mask_non_zeros.sum(axis=0)
                n_non_missing = np.add(n_non_missing, s)

            else:
                sums = X.sum(axis=axis)
                n_non_missing = np.diff(X.indptr)

            # Ignore the error, columns with a np.nan statistics_
            # are not an error at this point. These columns will
            # be removed in transform
            with np.errstate(all="ignore"):
                return np.ravel(sums) / np.ravel(n_non_missing)

        # Median + Most frequent
        else:
            # Remove the missing values, for each column
            columns_all = np.hsplit(X.data, X.indptr[1:-1])
            mask_missing_values = _get_mask(X.data, missing_values)
            mask_valids = np.hsplit(np.logical_not(mask_missing_values),
                                    X.indptr[1:-1])

            # astype necessary for bug in numpy.hsplit before v1.9
            columns = [col[mask.astype(bool, copy=False)]
                       for col, mask in zip(columns_all, mask_valids)]

            # Median
            if strategy == "median":
                median = np.empty(len(columns))
                for i, column in enumerate(columns):
                    median[i] = _get_median(column, n_zeros_axis[i])

                return median

            # Most frequent
            elif strategy == "most_frequent":
                most_frequent = np.empty(len(columns))

                for i, column in enumerate(columns):
                    most_frequent[i] = _most_frequent(column,
                                                      0,
                                                      n_zeros_axis[i])

                return most_frequent

    def _dense_fit(self, X, strategy, missing_values, axis):
        """Fit the transformer on dense data."""
        X = check_array(X, force_all_finite=False)
        mask = _get_mask(X, missing_values)
        masked_X = ma.masked_array(X, mask=mask)

        # Mean
        if strategy == "mean":
            mean_masked = np.ma.mean(masked_X, axis=axis)
            # Avoid the warning "Warning: converting a masked element to nan."
            mean = np.ma.getdata(mean_masked)
            mean[np.ma.getmask(mean_masked)] = np.nan

            return mean

        # Median
        elif strategy == "median":
            if tuple(int(v) for v in np.__version__.split('.')[:2]) < (1, 5):
                # In old versions of numpy, calling a median on an array
                # containing nans returns nan. This is different is
                # recent versions of numpy, which we want to mimic
                masked_X.mask = np.logical_or(masked_X.mask,
                                              np.isnan(X))
            median_masked = np.ma.median(masked_X, axis=axis)
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
            if axis == 0:
                X = X.transpose()
                mask = mask.transpose()

            most_frequent = np.empty(X.shape[0])

            for i, (row, row_mask) in enumerate(zip(X[:], mask[:])):
                row_mask = np.logical_not(row_mask).astype(np.bool)
                row = row[row_mask]
                most_frequent[i] = _most_frequent(row, np.nan, 0)

            return most_frequent

    def transform(self, X):
        """Impute all missing values in X.

        Parameters
        ----------
        X : {array-like, sparse matrix}, shape = [n_samples, n_features]
            The input data to complete.
        """
        if self.axis == 0:
            check_is_fitted(self, 'statistics_')
            X = check_array(X, accept_sparse='csc', dtype=FLOAT_DTYPES,
                            force_all_finite=False, copy=self.copy)
            statistics = self.statistics_
            if X.shape[1] != statistics.shape[0]:
                raise ValueError("X has %d features per sample, expected %d"
                                 % (X.shape[1], self.statistics_.shape[0]))

        # Since two different arrays can be provided in fit(X) and
        # transform(X), the imputation data need to be recomputed
        # when the imputation is done per sample
        else:
            X = check_array(X, accept_sparse='csr', dtype=FLOAT_DTYPES,
                            force_all_finite=False, copy=self.copy)

            if sparse.issparse(X):
                statistics = self._sparse_fit(X,
                                              self.strategy,
                                              self.missing_values,
                                              self.axis)

            else:
                statistics = self._dense_fit(X,
                                             self.strategy,
                                             self.missing_values,
                                             self.axis)

        # Delete the invalid rows/columns
        invalid_mask = np.isnan(statistics)
        valid_mask = np.logical_not(invalid_mask)
        valid_statistics = statistics[valid_mask]
        valid_statistics_indexes = np.where(valid_mask)[0]
        missing = np.arange(X.shape[not self.axis])[invalid_mask]

        if self.axis == 0 and invalid_mask.any():
            if self.verbose:
                warnings.warn("Deleting features without "
                              "observed values: %s" % missing)
            X = X[:, valid_statistics_indexes]
        elif self.axis == 1 and invalid_mask.any():
            raise ValueError("Some rows only contain "
                             "missing values: %s" % missing)

        # Do actual imputation
        if sparse.issparse(X) and self.missing_values != 0:
            mask = _get_mask(X.data, self.missing_values)
            indexes = np.repeat(np.arange(len(X.indptr) - 1, dtype=np.int),
                                np.diff(X.indptr))[mask]

            X.data[mask] = valid_statistics[indexes].astype(X.dtype,
                                                            copy=False)
        else:
            if sparse.issparse(X):
                X = X.toarray()

            mask = _get_mask(X, self.missing_values)
            n_missing = np.sum(mask, axis=self.axis)
            values = np.repeat(valid_statistics, n_missing)

            if self.axis == 0:
                coordinates = np.where(mask.transpose())[::-1]
            else:
                coordinates = mask

            X[coordinates] = values

        return X


class KNNImputer(BaseEstimator, TransformerMixin):
    """Imputation for completing missing values using k-Nearest Neighbors.

    Each sample's missing values are imputed from up to ``max_neighbors``
    nearest neighbors found in the training set. Each missing feature is then
    imputed as the average, either weighted or unweighted, of these neighbors
    who have a value for it. Where all neighbors have that feature value
    missing, the training set average for that feature is used for imputation.

    Parameters
    ----------
    missing_values : integer or "NaN", optional (default = "NaN")
        The placeholder for the missing values. All occurrences of
        `missing_values` will be imputed. For missing values encoded as
        ``np.nan``, use the string value "NaN".

    max_neighbors : int, optional (default = 5)
        Maximum number of neighboring samples to use for imputation. When any
        of the neighbors themselves have the feature value missing then the
        remaining neighbors, if any, that have the feature value available are
        used. But if none of the neighbors have the value available, the global
        feature mean (i.e., by default, the column mean) is used for
        imputation.

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
        The maximum percentage of columns (i.e. features) that can be missing
        before the sample is excluded from nearest neighbor imputation. It
        means that such rows will not be considered a potential donor in
        ``fit()``, and in ``transform()`` their missing feature values will be
        imputed to be the column mean for the entire dataset.

    col_max_missing : float, optional (default = 0.8)
        The maximum percentage of rows (or samples) that can be missing
        for a given feature beyond which an error is raised.

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
    >>> from sklearn.preprocessing.imputation import KNNImputer
    >>> nan = float("NaN")
    >>> X = [[1, 2, nan], [3, 4, 3], [nan, 6, 5], [8, 8, 7]]
    >>> imputer = KNNImputer(max_neighbors=2, weights="uniform")
    >>> imputer.fit_transform(X)
    array([[ 1. ,  2. ,  4. ],
           [ 3. ,  4. ,  3. ],
           [ 5.5,  6. ,  5. ],
           [ 8. ,  8. ,  7. ]])
    """

    def __init__(self, missing_values="NaN", max_neighbors=5,
                 weights="uniform", metric="masked_euclidean",
                 row_max_missing=0.5, col_max_missing=0.8,
                 use_complete=False, copy=True):

        self.missing_values = missing_values
        self.max_neighbors = max_neighbors
        self.weights = weights
        self.metric = metric
        self.row_max_missing = row_max_missing
        self.col_max_missing = col_max_missing
        self.use_complete = use_complete
        self.copy = copy

    def fit(self, X, y=None):
        """Fit the imputer on X.

        Parameters
        ----------
        X : {array-like}, shape (n_samples, n_features)
            Input data, where ``n_samples`` is the number of samples and
            ``n_features`` is the number of features.

        Returns
        -------
        self : object
            Returns self.
        """
        # Imports here to avoid circular import
        from ..neighbors import NearestNeighbors
        from ..neighbors.base import _check_weights

        # Check parameters
        force_all_finite = False if self.missing_values in ["NaN",
                                                            np.nan] else True
        X = check_array(X, accept_sparse=False, dtype=np.float64,
                        force_all_finite=force_all_finite, copy=self.copy)
        self.weights = _check_weights(self.weights)

        # Check for +/- inf
        if np.any(np.isinf(X)):
            raise ValueError("+/- inf values are not allowed.")

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
        if X.shape[0] < self.max_neighbors:
            raise ValueError("There are only %d samples, but max_neighbors=%d."
                             % (X.shape[0], self.max_neighbors))

        # Instantiate NN object, get column means, and store in statistics_
        neigh = NearestNeighbors(n_neighbors=self.max_neighbors,
                                 metric=self.metric,
                                 metric_params={"missing_values":
                                                self.missing_values})
        self._fitted_neighbors = neigh.fit(X)
        self.statistics_ = X_col_means

        return self

    def _get_weight_matrix(self, fitted_X, mask, adjusted_max_neighbors,
                           receiver_row_index, row_repeats,
                           knn_row_index, knn_distances):
        """Get the weight matrix for the donors"""

        # Import(s) here to avoid circular import
        from ..neighbors.base import _get_weights

        # If different X in transform, get a new mask
        if self.max_neighbors == adjusted_max_neighbors:
            nbors_mask = _get_mask(fitted_X[knn_row_index],
                                   value_to_mask=self.missing_values)
        else:
            nbors_mask = mask[knn_row_index]

        # Anti-mask tells us what is NOT missing
        nbors_anti_mask = ~nbors_mask
        receiver_anti_mask = ~mask[receiver_row_index]

        # Sum anti-masks to see if both donor & receiver are missing
        # A zero value indicates that a feature is missing in both
        # Sum over all cols to locate degenerate donors
        anti_masks_combined = receiver_anti_mask + nbors_anti_mask
        anti_masks_combined = anti_masks_combined.sum(axis=-1)
        degenerate_nbors = anti_masks_combined < mask.shape[1]
        knn_distances[degenerate_nbors] = np.inf

        # Retreive and, if applicable, transform weight matrix
        weight_matrix = _get_weights(knn_distances, self.weights)
        if weight_matrix is not None:
            weight_matrix = weight_matrix[:, np.newaxis, :]
            weight_matrix = np.repeat(weight_matrix,
                                      row_repeats, axis=0).ravel()
            weight_matrix = weight_matrix.reshape(
                (-1, adjusted_max_neighbors))
        return weight_matrix

    def _transform(self, X, adjusted_max_neighbors):
        """Impute all missing values in X.

        Parameters
        ----------
        X : {array-like}, shape = [n_samples, n_features]
            The input data to complete.

        adjusted_max_neighbors : int
            Depending on the calling method, the default value must
            either be equal to max_neighbors or max_neighbors + 1.
            If the calling method is transform(), then its value needs to be
            equal to max_neighbors and if calling method is fit_transform()
            then its value must be equal to max_neighbors + 1.
        """

        check_is_fitted(self, 'statistics_')
        force_all_finite = False if self.missing_values in ["NaN",
                                                            np.nan] else True
        X = check_array(X, accept_sparse=False, dtype=FLOAT_DTYPES,
                        force_all_finite=force_all_finite, copy=self.copy)
        # Check for +/- inf
        if np.any(np.isinf(X)):
            raise ValueError("+/- inf values are not allowed in data to be "
                             "transformed.")

        # Get fitted data and ensure correct dimension
        fitted_X = self._fitted_neighbors._fit_X
        if X.shape[1] != fitted_X.shape[1]:
            raise ValueError("Incompatible dimension between the fitted "
                             "dataset and the one to be transformed.")
        mask = _get_mask(X, self.missing_values)
        n_rows_X, n_cols_X = X.shape
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
            if self.use_complete:
                # Initializations
                # Mask for fitted_X
                mask_fx = _get_mask(fitted_X, np.nan)

                # Locate unique patterns, but first a numpy version check
                np_version = np.__version__.split(".")

                # Different behavior for np.unique in Numpy < 1.13.x
                if float('.'.join(np_version[:2])) < 1.13:
                    patterns, row_pat_idx = _unique(
                        mask, return_inverse=True, axis=0)
                else:
                    patterns, row_pat_idx = np.unique(
                        mask, return_inverse=True, axis=0)

                # Get row idx for receivers (missing)
                receiver_row_missing_index, _ = np.where(mask)

                # For every pattern, index receivers and potential donors
                for p in range(len(patterns)):
                    if not np.any(patterns[p]):
                        continue
                    # receivers are those with pattern 'p'
                    row_has_missing_pat = (row_pat_idx == p)

                    # Donors have features missing in receivers available
                    # The bitwise-AND captures if something is missing in both
                    donor_row_idx = ~np.any(patterns[p] & mask_fx, axis=1)

                    # Change donor set to ones not missing relevant features
                    self._fitted_neighbors._fit_X = fitted_X[donor_row_idx]
                    if len(self._fitted_neighbors._fit_X) < self.max_neighbors:
                        err_msg = "Insufficient neighbors with feature values."
                        raise ValueError(err_msg)

                    # Row index of receivers & identify potential donors
                    receiver_row_index = np.where(
                        row_has_missing_pat)[0].reshape((-1, 1))
                    neighbors = self._fitted_neighbors.kneighbors(
                        X[row_has_missing_pat, :],
                        n_neighbors=self.max_neighbors)

                    # Get row index, distance, and weights of donors
                    knn_distances, knn_row_index = neighbors
                    row_repeats = row_total_missing[row_has_missing_pat]

                    # Weighting: Set self/degenerate donor(s) distance to inf
                    weight_matrix = None
                    if self.weights in ["distance"] or callable(self.weights):
                        weight_matrix = self._get_weight_matrix(
                            self._fitted_neighbors._fit_X,
                            mask,
                            self.max_neighbors,
                            receiver_row_index,
                            row_repeats,
                            knn_row_index,
                            knn_distances
                        )

                    # Repeat each set donor indices by
                    # missing count in the corresponding recipient row
                    knn_row_index_repeat = np.repeat(
                        knn_row_index, row_repeats, axis=0).ravel()

                    # Get repeated column index of donors
                    _, receiver_col_missing_index = \
                        np.where(mask[row_has_missing_pat])
                    knn_col_index_repeat = np.repeat(
                        receiver_col_missing_index,
                        self.max_neighbors)

                    # Retrieve donor cells and calculate kNN score
                    donors = self._fitted_neighbors._fit_X[
                        knn_row_index_repeat, knn_col_index_repeat].reshape(
                        (-1, self.max_neighbors))
                    donors_mask = _get_mask(donors, self.missing_values)
                    donors = np.ma.array(donors, mask=donors_mask)

                    # Final imputation
                    imputed = np.ma.average(donors, axis=1,
                                            weights=weight_matrix)
                    X[np.where(row_has_missing_pat)[0],
                      receiver_col_missing_index] = imputed.data

                    # Recover original dataset
                    self._fitted_neighbors._fit_X = fitted_X
            else:
                # Row index of receivers & identify potential donors
                receiver_row_index = np.where(
                    row_has_missing)[0].reshape((-1, 1))
                neighbors = self._fitted_neighbors.kneighbors(
                    X[row_has_missing, :], n_neighbors=adjusted_max_neighbors)

                # Get row index, distance, and weights of donors
                knn_distances, knn_row_index = neighbors
                row_repeats = row_total_missing[row_total_missing != 0]

                # Weighting: Set self and degenerate donor(s) distance to inf
                weight_matrix = None
                if self.weights in ["distance"] or callable(self.weights):
                    weight_matrix = self._get_weight_matrix(
                        fitted_X,
                        mask,
                        adjusted_max_neighbors,
                        receiver_row_index,
                        row_repeats,
                        knn_row_index,
                        knn_distances
                    )

                # Repeat each set donor indices by
                # missing count in the corresponding recipient row
                knn_row_index_repeat = np.repeat(
                    knn_row_index, row_repeats, axis=0).ravel()

                # Get repeated column index of donors
                receiver_row_missing_index, receiver_col_missing_index = \
                    np.where(mask)
                knn_col_index_repeat = np.repeat(receiver_col_missing_index,
                                                 adjusted_max_neighbors)

                # Retrieve donor cells and calculate kNN score
                donors = fitted_X[
                    knn_row_index_repeat, knn_col_index_repeat].reshape(
                    (-1, adjusted_max_neighbors))
                donors_mask = _get_mask(donors, self.missing_values)
                donors = np.ma.array(donors, mask=donors_mask)

                # Warning if donor count < max_neighbors
                if np.any(donors_mask.sum(axis=1) < self.max_neighbors):
                    warnings.warn("One or more donor(s) have the relevant "
                                  "feature value missing.")

                # Final imputation
                imputed = np.ma.average(donors, axis=1, weights=weight_matrix)
                X[mask] = imputed.data
                unimputed_index = np.where(donors_mask.all(axis=1))
                if len(unimputed_index[0]) > 0:
                    unimputed_rows = receiver_row_missing_index[
                        unimputed_index]
                    unimputed_cols = receiver_col_missing_index[
                        unimputed_index]
                    X[unimputed_rows, unimputed_cols] = np.take(
                        self.statistics_, unimputed_cols)

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

        This method should *only* be used if the data to be fitted is the
        same as the data to be transformed.

        Parameters
        ----------
        X : {array-like}, shape (n_samples, n_features)
            Input data, where ``n_samples`` is the number of samples and
            ``n_features`` is the number of features.

        Returns
        -------
        X : {array-like}, shape (n_samples, n_features)
            Returns imputed dataset.
        """
        return self.fit(X)._transform(
            X, adjusted_max_neighbors=self.max_neighbors + 1)

    def transform(self, X):
        """Impute all missing values in X.

        This method should *only* be used if the data to be fitted is different
        from the data to be transformed.

        WARNING: If the same dataset is passed in fit() and transform(),
        one of the returned "neighbors" maybe the sample itself. If you will be
        passing the same dataset, use fit_transform() to avoid this behavior.

        Parameters
        ----------
        X : {array-like}, shape = [n_samples, n_features]
            The input data to complete.

        Returns
        -------
        X : {array-like}, shape (n_samples, n_features)
            Returns imputed dataset.
        """
        check_is_fitted(self, 'statistics_')
        return self._transform(X, adjusted_max_neighbors=self.max_neighbors)
