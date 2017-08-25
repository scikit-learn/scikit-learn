# Authors: Nicolas Tresegnie <nicolas.tresegnie@gmail.com>
# License: BSD 3 clause

import warnings

import numpy as np
import numpy.ma as ma
from scipy import sparse
from scipy import stats

from ..base import BaseEstimator, TransformerMixin
from ..utils import check_array
from ..utils import deprecated
from ..utils.sparsefuncs import _get_median
from ..utils.validation import check_is_fitted
from ..utils.validation import FLOAT_DTYPES

from ..externals import six

zip = six.moves.zip
map = six.moves.map

__all__ = [
    'Imputer', 'MissingIndicator'
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


@deprecated("Imputer was deprecated in version 0.20 and will be "
            "removed in 0.22. Import impute.SimpleImputer from "
            "sklearn instead.")
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


class MissingIndicator(BaseEstimator, TransformerMixin):
    """Binary indicators for missing values.

    Parameters
    ----------
    missing_values : integer or np.NaN, optional (default=np.NaN)
        The placeholder for the missing values. All occurrences of
        ``missing_values`` will be represented as ones.

    features : {'auto' (default), 'all', array-like of int}
        If "all", mask will represent all features.
        If "auto", mask will only represent features with missing values
        during fit time.
        If mask/indices, mask will only represent features in the
        indices or mask.

    sparse : boolean or "auto", optional (default="auto")
        If True, the transformed ``X`` will be a sparse matrix.
        If False, the transformed ``X`` will be a numpy array.
        If "auto", the transformed ``X`` will be of same type as input.

    error_on_new : boolean, optional (default=True)
        If True, transform will raise an error when there are features with
        missing values in transform but have no missing values in fit
        This is applicable only when ``features`` = "auto"

    Example
    -------
    >>> from sklearn.preprocessing import MissingIndicator
    >>> import numpy as np
    >>> X1 = np.array([
    ...       ["NaN",  1,  3],
    ...       [ 4,  0, "NaN"],
    ...       [ 8,  1,  0]
    ... ])
    >>> X2 = np.array([
    ...       [ 5, 1, "NaN"],
    ...       ["NaN",  2,  3],
    ...       [ 2,  4,  0]
    ... ])
    >>> indicator = MissingIndicator()
    >>> indicator.fit(X1)
    MissingIndicator(error_on_new=True, features='auto', missing_values='NaN',
             sparse='auto')
    >>> X2_tr = indicator.transform(X2)
    >>> X2_tr
    array([[0, 1],
           [1, 0],
           [0, 0]])

    Attributes
    ----------
    feat_with_missing_  : array of shape(n_missing_features,)
        The features with missing values.

    n_features_ : int
        The number of features in the input.
    """

    def __init__(self, missing_values="NaN", features="auto", sparse="auto",
                 error_on_new=True):
        self.missing_values = missing_values
        self.features = features
        self.sparse = sparse
        self.error_on_new = error_on_new

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
        self._validate_params()

        X = check_array(X, accept_sparse=('csc', 'csr'), dtype=np.float64, force_all_finite=False)
        self.n_features_ = X.shape[1]

        if isinstance(self.features, six.string_types):
            if self.features == "auto":
                _, self.feat_with_missing_ = self._get_missing_features_info(X)
            else:  # self.features == "all"
                self.feat_with_missing_ = np.arange(self.n_features_)
        else:
            self.feat_with_missing_ = self.features

        return self

    def transform(self, X):
        """Generate missing values indicator for X.

        Parameters
        ----------
        X : {array-like, sparse matrix}, shape = [n_samples, n_features]
            The input data to complete.

        Returns
        -------
        Xt : array or sparse matrix, shape = [n_samples, n_features]
             The missing indicator for input data

        """
        check_is_fitted(self, "feat_with_missing_", "n_features_")

        X = check_array(X, accept_sparse=('csc', 'csr'), dtype=np.float64, force_all_finite=False)
        if X.shape[1] != self.n_features_:
            raise ValueError("X has a different number of features "
                             "than during fitting.")

        imputer_mask, feat_with_missing = self._get_missing_features_info(X)

        if isinstance(self.features, six.string_types):
            if self.features == "auto":
                print feat_with_missing
                print self.feat_with_missing_
                features = np.setdiff1d(feat_with_missing,
                                        self.feat_with_missing_)
                if self.error_on_new and features.size > 0:
                    raise Exception("The features %s have missing values "
                                    "in transform but have no missing values "
                                    "in fit" % features)

        if not (isinstance(self.features, six.string_types) and self.features == "all"):
            # no need to slice when all features have missing values
            imputer_mask = imputer_mask[:, self.feat_with_missing_]

        return imputer_mask

    def fit_transform(self, X, y=None):
        """Generate missing values indicator for X.

        Parameters
        ----------
        X : {array-like, sparse matrix}, shape = [n_samples, n_features]
            The input data to complete.

        Returns
        -------
        Xt : array or sparse matrix, shape = [n_samples, n_features]
             The missing indicator for input data

        """
        self.fit(X, y)
        return self.transform(X)

    def _validate_params(self):
        if (isinstance(self.features, six.string_types) and
                self.features not in ["auto", "all"]):
            raise ValueError("Can only use these options: 'auto', 'all'"
                             " got {0}".format(self.features))
        elif not isinstance(self.features, six.string_types):
            self.features = check_array(self.features, ensure_2d=False)
            if (isinstance(self.features, np.ndarray) and
                    self.features.dtype.kind != 'i'):
                raise ValueError("Features should be an array of integers")

        if not ((isinstance(self.sparse, six.string_types) and
                self.sparse == "auto") or isinstance(self.sparse, bool)):
            raise ValueError("sparse can only boolean or 'auto'"
                             " got {0}".format(self.sparse))

    def _get_missing_features_info(self, X):
        if sparse.issparse(X) and self.missing_values != 0:
            #  sparse matrix and missing values is not zero
            imputer_mask = _get_mask(X.data, self.missing_values)
            imputer_mask = X.__class__((imputer_mask, X.indices.copy(),
                                        X.indptr.copy()), shape=X.shape,
                                       dtype=X.dtype)
            feat_with_missing = imputer_mask.sum(axis=0).nonzero()[1]
        else:
            #  sparse with zero as missing value and dense matrix
            if sparse.issparse(X):
                X = X.toarray()
            imputer_mask = _get_mask(X, self.missing_values)
            #  convert boolean mask to binary mask
            imputer_mask = imputer_mask.astype(int, copy=False)
            feat_with_missing = imputer_mask.sum(axis=0).nonzero()[0]

        if ((self.sparse == 'auto' and sparse.issparse(imputer_mask)) or
                self.sparse is True):
            imputer_mask = sparse.csc_matrix(imputer_mask)
        elif self.sparse is False and sparse.issparse(imputer_mask):
            imputer_mask = imputer_mask.toarray()

        print 'before'
        print feat_with_missing
        feat_with_missing = feat_with_missing.ravel()
        print 'after'
        print feat_with_missing
        return imputer_mask, feat_with_missing
