# Authors: Nicolas Tresegnie <nicolas.tresegnie@gmail.com>
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
# from ..neighbors.base import _get_weights, _check_weights

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
        self : object
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
    """Imputation for completing missing values using Nearest Neighbors.

    Broadly speaking, the imputation is performed using either
    the weighted or the unweighted mean of the desired number of neighbors.

    Parameters
    ----------
    missing_values : integer or "NaN", optional (default="NaN")
        The placeholder for the missing values. All occurrences of
        `missing_values` will be imputed. For missing values encoded as np.nan,
        use the string value "NaN".

    n_neighbors : int, optional (default = 5)
        Maximum number of neighboring samples to use for imputation. When
        any of the neighbors themselves have the feature value missing then
        the remaining n_neighbors-1 neighbors are used and, if need be,
        the process repeats until a single neighbor remains. If all
        the neighbors have the feature value missing, then the overall feature
        mean is used for imputation.

    weights : str or callable
        weight function used in prediction.  Possible values:

        - 'uniform' : uniform weights.  All points in each neighborhood
          are weighted equally.
        - 'distance' : weight points by the inverse of their distance.
          in this case, closer neighbors of a query point will have a
          greater influence than neighbors which are further away.
        - [callable] : a user-defined function which accepts an
          array of distances, and returns an array of the same shape
          containing the weights.

        Uniform weights are used by default.

    metric : string or callable, optional (default = 'masked_euclidean')
        metric to use for distance computation.

    row_max_missing : float, optional (default = 0.5)
        The maximum percentage of columns (i.e. features) that can be missing
        before the sample is excluded from nearest neighbor imputation. It
        means that such rows will not be considered a potential donor in fit()
        and in transform() their missing feature values will be imputed to be
        the column mean for the entire dataset.

    col_max_missing : float, optional (default = 0.8)
        The maximum percentage of rows (or samples) that can be missing
        for a given feature beyond which an error is raised.

    copy : boolean, optional (default=True)
        If True, a copy of X will be created. If False, imputation will
        be done in-place whenever possible. Note that, if metric is
        "masked_euclidean" and copy=False then missing_values in the
        input matrix X will be overwritten with zeros.

    Attributes
    ----------
    statistics_ : {tuple}
        A tuple whose first element is the input dataset used to fit the
        KNNImputer object and the second element is the column means of that
        dataset using observed (i.e. non-missing) values.

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
    >>> imputer = KNNImputer(n_neighbors=2, weights="uniform")
    >>> imputer.fit_transform(X)
    array([[ 1. ,  2. ,  4. ],
           [ 3. ,  4. ,  3. ],
           [ 5.5,  6. ,  5. ],
           [ 8. ,  8. ,  7. ]])
    """

    def __init__(self, missing_values="NaN", n_neighbors=5,
                 weights="uniform", metric="masked_euclidean",
                 row_max_missing=0.5, col_max_missing=0.8, copy=True):

        self.missing_values = missing_values
        self.n_neighbors = n_neighbors
        self.weights = weights
        self.metric = metric
        self.row_max_missing = row_max_missing
        self.col_max_missing = col_max_missing
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
        if (np.any(np.isinf(X))):
            raise ValueError("+/- inf values are not allowed.")

        # Check if % missing in any column > col_max_missing
        mask = _get_mask(X, self.missing_values)
        if np.any(mask.sum(axis=0) > (X.shape[0] * self.col_max_missing)):
            raise ValueError("The following columns have, "
                             "more than {0}% missing values: {1}"
                             .format(self.col_max_missing*100, np.where(
                              mask.sum(axis=0) > (X.shape[0] *
                                                  self.col_max_missing))))
        X_col_means = np.ma.array(X, mask=mask).mean(axis=0).data

        # Check if % missing in any row > col_max_missing
        bad_rows = mask.sum(axis=1) > (mask.shape[1] * self.row_max_missing)
        if np.any(bad_rows):
            warnings.warn(
                "The following rows have more than {0}% missing values and "
                "are not included as donor neighbors: {1}"
                    .format(self.row_max_missing*100, np.where(bad_rows)))

            # Remove rows that have more than row_max_missing % missing
            X = X[~bad_rows, :]

        # Check if sufficient neighboring samples available
        if X.shape[0] < self.n_neighbors:
            raise ValueError("There are only %d samples, "
                             "but n_neighbors=%d."
                             % (X.shape[0], self.n_neighbors))

        # Instantiate NN object, get column means, and store in statistics_
        neigh = NearestNeighbors(n_neighbors=self.n_neighbors,
                                 metric=self.metric)
        self._fitted_neighbors = neigh.fit(X)
        self.statistics_ = X_col_means

        return self

    def _transform(self, X, n_neighbors_new):
        """Impute all missing values in X.

        Parameters
        ----------
        X : {array-like}, shape = [n_samples, n_features]
            The input data to complete.

        n_neighbors_new : int
            Indicates whether to pass n_neighbors or n_neighbors+1 to
            _tranform().
            Calling transform() automatically sets this to self.n_neighbors
            while fit_transform() sets it to self.n_neighbors + 1.
        """
        # Import(s) here to avoid circular import
        from ..neighbors.base import _get_weights

        check_is_fitted(self, 'statistics_')
        force_all_finite = False if self.missing_values in ["NaN",
                                                            np.nan] else True
        X = check_array(X, accept_sparse=False, dtype=FLOAT_DTYPES,
                        force_all_finite=force_all_finite, copy=self.copy)
        # Check for +/- inf
        if (np.any(np.isinf(X))):
            raise ValueError("+/- inf values are not allowed.")

        # Get fitted data and ensure correct dimension
        fitted_X = self._fitted_neighbors._fit_X
        if X.shape[1] != fitted_X.shape[1]:
            raise ValueError("Incompatible dimension between the fitted "
                             "dataset and the one to be transformed.")
        mask = _get_mask(X, self.missing_values)
        n_rows_X, n_cols_X = X.shape
        row_total_missing = mask.sum(axis=1)
        if not np.any(row_total_missing > 0):
            return X
        # row_has_missing = row_total_missing.astype(np.bool)

        # Check for excessive missingness in rows
        bad_rows = row_total_missing > (mask.shape[1] * self.row_max_missing)
        if np.any(bad_rows):
            warnings.warn(
                "The following rows have more than {0}% missing values and "
                "are imputed with column means: {1}"
                    .format(self.row_max_missing*100, np.where(bad_rows)))
            X_bad = X[bad_rows, :]
            X = X[~bad_rows, :]
            mask = _get_mask(X, self.missing_values)
            row_total_missing = mask.sum(axis=1)
        row_has_missing = row_total_missing.astype(np.bool)
        # Check if the X in fit() and transform() are the same
        # if X is fitted_X or X.base is fitted_X.base:
        #     neighbors = self._fitted_neighbors.kneighbors(
        #         n_neighbors=self.n_neighbors)
        # else:
        #     neighbors = self._fitted_neighbors.kneighbors(
        #         X[row_has_missing, :], n_neighbors=self.n_neighbors)
        if np.any(row_has_missing):
            neighbors = self._fitted_neighbors.kneighbors(
                X[row_has_missing, :], n_neighbors=n_neighbors_new)

            # Get row index, distance, and weights of donors
            knn_distances, knn_row_index = neighbors
            if n_neighbors_new > self.n_neighbors:
                knn_distances = knn_distances[:, 1:]
                knn_row_index = knn_row_index[:, 1:]
            weights = _get_weights(knn_distances, self.weights)

            knn_row_index = np.vsplit(knn_row_index, knn_row_index.shape[0])
            row_repeats = row_total_missing[row_total_missing != 0]
            knn_row_index = np.repeat(
                knn_row_index, row_repeats, axis=0).ravel()

            # Get column index of donors
            # NOTE: Following assumes columns in X and _fit_X are in the same
            # order
            row_missing_index, col_missing_index = np.where(mask)
            knn_col_index = np.repeat(col_missing_index, self.n_neighbors)

            # Calculate kNN score and impute
            donors = fitted_X[
                (knn_row_index, knn_col_index)].reshape((-1, self.n_neighbors))
            donors_mask = _get_mask(donors, self.missing_values)
            donors = np.ma.array(
                donors, mask=donors_mask)
            imputed = np.ma.average(donors, axis=1, weights=weights)
            X[mask] = imputed.data
            unimputed_index = np.where(
                donors_mask.sum(axis=1) == self.n_neighbors)
            # imputed_mask = _get_mask(imputed.data, self.missing_values)
            if len(unimputed_index[0]) > 0:
                # unimputed_loc = np.where(imputed_mask)
                # unimputed_rows, unimputed_cols = np.where(mask)
                unimputed_rows = row_missing_index[unimputed_index]
                unimputed_cols = col_missing_index[unimputed_index]
                X[(unimputed_rows, unimputed_cols)] = np.take(self.statistics_,
                                                              unimputed_cols)

        # Merge bad rows to X and mean impute any leftover missing
        if np.any(bad_rows):
            bad_missing_index = np.where(_get_mask(X_bad, self.missing_values))
            X_bad[bad_missing_index] = np.take(self.statistics_,
                                               bad_missing_index[1])
            X_merged = np.empty((n_rows_X, n_cols_X))
            X_merged[bad_rows, :] = X_bad
            X_merged[~bad_rows, :] = X
            X = X_merged

        # Impute bad_rows and leftover missing with column means
        # mask_after_knn = _get_mask(X[row_has_missing, :],
        #                            self.missing_values)
        # if np.any(mask_after_knn):
        #     missing_index = np.where(mask_after_knn)
        #     X[row_has_missing, :][missing_index] = np.take(self.statistics_,
        #                                                    missing_index[1])

        return X

    def fit_transform(self, X, y=None, **fit_params):
        """Impute all missing values in X.

        Parameters
        ----------
        X : {array-like}, shape = [n_samples, n_features]
            The input data to complete.
        """
        return self.fit(X)._transform(X, n_neighbors_new=self.n_neighbors + 1)

    def transform(self, X):
        """Impute all missing values in X.

        Parameters
        ----------
        X : {array-like}, shape = [n_samples, n_features]
            The input data to complete.
        """
        check_is_fitted(self, 'statistics_')
        return self._transform(X, n_neighbors_new=self.n_neighbors)
