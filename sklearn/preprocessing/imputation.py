# Authors: Nicolas Tresegnie <nicolas.tresegnie@gmail.com>
# License: BSD 3 clause

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

def _get_knn(X, value_to_mask="NaN", n_neighbors=10, Y=None):
    """Returns the k(=n_neighbors) nearest neighbors of vectors in a
     given matrix in euclidean space. If two matrices are passed,
     then k-Nearest Neighbors of vectors in X in the matrix Y is returned."""

    # Setup missing mask
    mask_X = _get_mask(X, value_to_mask)

    if Y is None:
        # Setup the anti-mask and change missing to zero
        mask_X = _get_mask(X, value_to_mask)
        XT = np.transpose(X)
        N = (~mask_X) * 1
        NT = np.transpose(N)
        X[mask_X] = 0

        # Matrix formula to calculate pair-wise distance between all vectors in a
        # matrix with missing values. It zero-weights coordinates with missing value
        # in either vector in the pair and up-weights the remaining coordinates.
        # Matrix formula derived by: Shreya Bhattarai <shreya.bhattarai@gmail.com>

        """
        Store np.dot(N, (XT * XT)) and add its transpose rather than 
        redoing a matrix product
        dist = np.sqrt((X.shape[1] * 1 / ((np.dot(N, NT)))) * (
            np.dot(N, (XT * XT)) - 2 * (np.dot(X, XT)) +
            np.dot((X * X), NT)))        

        N_dot_XT2 = np.dot(N, (XT * XT))
        N_dot_XT2_T = np.transpose(N_dot_XT2)
        """

        N_dot_XT2 = np.dot(N, (XT * XT))
        N_dot_XT2_T = np.transpose(N_dot_XT2)

        dist = np.sqrt((X.shape[1] * 1 / ((np.dot(N, NT)))) * (
            N_dot_XT2 - 2 * (np.dot(X, XT)) +
            N_dot_XT2_T))

        # Set distance with self to np.inf
        np.fill_diagonal(dist, np.inf)

    else:
        # ValueError if X and Y have incompatible dimensions
        if X.shape[1] != Y.shape[1]:
            raise ValueError("The search dimension of the matrices "
                             "are not equal: [{0}] versus [{1}]".
                             format(X.shape[1], Y.shape[1]))

        mask_Y = _get_mask(Y, value_to_mask)
        NY = (~mask_Y) * 1
        YT = np.transpose(Y)
        mask_YT = _get_mask(YT, value_to_mask)
        NYT = np.transpose(NY)
        YT[mask_YT] = 0

        NX = (~mask_X) * 1
        X[mask_X] = 0

        # Matrix formula to calculate pair-wise distance between all vectors in a
        # matrix X to vectors in matrix Y. It handles missing values the same way
        # as for a single matrix.
        # Matrix formula derived by: Shreya Bhattarai <shreya.bhattarai@gmail.com>

        dist = np.sqrt((X.shape[1] * 1 / ((np.dot(NX, NYT)))) *
                       (np.dot((X * X), NYT) - 2 * (np.dot(X, YT)) + np.dot(NX, (YT * YT))))

    # Ensure enough candidate neighbors are available
    n_candidates = dist.shape[1] if Y is not None else dist.shape[1] - 1
    if n_candidates < n_neighbors:
        raise ValueError("There are only %d candidate neighbors, "
                         "but n_neighbors=%d."
                         % (dist.shape[1] - 1, n_neighbors))

    # Missing locations and counts
    row_missing_sum_X = mask_X.sum(axis=1)
    # is_row_missing_X = np.any(mask_X, axis=1)
    # is_col_missing_X = np.any(mask_X, axis=0)
    col_missing_index_X = np.where(mask_X)[1]

    # Arg-partition (quasi-argsort) of n_neighbors and retrieve them
    nbors_index = np.argpartition(dist, n_neighbors - 1, axis=1)
    knn_row_index = nbors_index[:, :n_neighbors]
    knn_row_index = np.vsplit(knn_row_index, knn_row_index.shape[0])
    knn_row_index = np.repeat(knn_row_index, row_missing_sum_X, axis=0)
    knn_row_index = knn_row_index.ravel()
    # This assumes columns in X and Y are in the same order; maybe change this?
    knn_col_index = np.repeat(col_missing_index_X, n_neighbors)
    return knn_row_index, knn_col_index

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
                 axis=0, verbose=0, copy=True, n_neighbors=10,
                 row_max_missing=0.5, col_max_missing=0.8):
        self.missing_values = missing_values
        self.strategy = strategy
        self.axis = axis
        self.verbose = verbose
        self.copy = copy
        self.n_neighbors = n_neighbors
        self.row_max_missing = row_max_missing
        self.col_max_missing = col_max_missing

    def fit(self, X, y=None):
        """Fit the imputer on X.

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
        # Check parameters
        allowed_strategies = ["mean", "median", "most_frequent", "knn"]
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

            # KNN
            elif strategy == "knn":
                raise ValueError("strategy='knn' does not support sparse "
                                 "matrix input yet.")

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

        # KNN
        elif strategy == "knn":
            if self.copy:
                X = np.copy(X)

            if axis == 1:
                X = X.transpose()
                mask = mask.transpose()

            #Get dimensions and missing count
            n_rows, n_cols = X.shape
            row_missing_sum = mask.sum(axis=1)

            #ValueError if % missing in any column > self.col_max_missing
            if np.any(mask.sum(axis=0) > (X.shape[0] * self.col_max_missing)):
                raise ValueError("The following axis position(s) have, "
                                 "more than {0}% missing values: {1}"
                                 .format(self.col_max_missing*100,np.where(mask.sum(axis=0) >
                                                  (X.shape[0] * self.col_max_missing))))

            if X.shape[0] < self.n_neighbors:
                raise ValueError("There are only %d samples, "
                                 "but n_neighbors=%d."
                                 % (X.shape[0], self.n_neighbors))

            #Fit to data

            # Check for excessive missingness in rows
            bad_rows = row_missing_sum > (mask.shape[1] * self.row_max_missing)
            X_bad = X[bad_rows, :]

            if np.any(bad_rows):
                X = X[~bad_rows, :]
                mask = _get_mask(X, missing_values)

            #Get the k nearest neighbors and impute
            if hasattr(self, 'statistics_'):
                Y = self.statistics_.data
                knnrows_index, knncols_index = _get_knn(X,
                                                        n_neighbors=self.n_neighbors, Y=Y)
                X[mask] = np.nan
                imputed = np.nanmean((Y[(knnrows_index, knncols_index)]).
                                     reshape((-1, self.n_neighbors)), axis=1)
            else:
                knnrows_index, knncols_index = _get_knn(X,
                                                        n_neighbors=self.n_neighbors)
                X[mask] = np.nan
                imputed = np.nanmean((X[(knnrows_index, knncols_index)]).
                                     reshape((-1, self.n_neighbors)), axis=1)
            X[mask] = imputed

            #Merge bad rows to X and mean impute any leftover missing
            if np.any(bad_rows):
                X_merged = np.empty((n_rows, n_cols))
                X_merged[bad_rows, :] = X_bad
                X_merged[~bad_rows, :] = X
                X = X_merged

            #Impute bad_rows and leftover missing with column means
            mask_after_knn = _get_mask(X, self.missing_values)
            if np.any(mask_after_knn):
                missing_index = np.where(mask_after_knn)
                X_col_means = masked_X.mean(axis=0).data
                X[missing_index] = np.take(X_col_means, missing_index[1])

            # Transpose back
            if axis == 1:
                X = X.transpose()

            #The mask is used to compare this imputed matrix with
            #input matrix in transform(), so return X as a masked array.
            X = np.ma.array(X,mask=masked_X.mask)
            return X

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
            #> Added knn exception below
            if self.strategy != "knn" and X.shape[1] != statistics.shape[0]:
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

        if self.strategy != "knn":
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

            if self.strategy == 'knn':
                if self.axis == 1:
                    X = X.transpose()
                    mask = mask.transpose()
                    statistics = statistics.transpose()

                #Check if the masks and the unmasked values are equal
                mask_fitted = statistics.mask
                masked_X = np.ma.array(X, mask=mask)
                if np.array_equal(mask, mask_fitted)\
                        and np.ma.allequal(masked_X, statistics):
                        X = statistics.data
                else:
                    X = self._dense_fit(X,
                                                       self.strategy,
                                                       self.missing_values,
                                                       self.axis).data

                if self.axis == 1:
                    X = X.transpose()

            else:
                values = np.repeat(valid_statistics, n_missing)

                if self.axis == 0:
                    coordinates = np.where(mask.transpose())[::-1]
                else:
                    coordinates = mask

                X[coordinates] = values

        return X
