"""Transformers for missing value imputation"""
# Authors: Nicolas Tresegnie <nicolas.tresegnie@gmail.com>
#          Sergey Feldman <sergeyfeldman@gmail.com>
# License: BSD 3 clause

from __future__ import division

import warnings
import numbers
import time

import numpy as np
import numpy.ma as ma
from scipy import sparse
from scipy import stats
from collections import namedtuple

from .base import BaseEstimator, TransformerMixin
from .base import clone
from .preprocessing import normalize
from .utils import check_array, check_random_state, safe_indexing
from .utils.sparsefuncs import _get_median
from .utils.validation import check_is_fitted
from .utils.validation import FLOAT_DTYPES
from .utils.fixes import _object_dtype_isnan
from .utils import is_scalar_nan

from .externals import six

zip = six.moves.zip
map = six.moves.map

ImputerTriplet = namedtuple('ImputerTriplet', ['feat_idx',
                                               'neighbor_feat_idx',
                                               'predictor'])

__all__ = [
    'MissingIndicator',
    'SimpleImputer',
    'ChainedImputer',
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


class ChainedImputer(BaseEstimator, TransformerMixin):
    """Chained imputer transformer to impute missing values.

    Basic implementation of chained imputer from MICE (Multivariate
    Imputations by Chained Equations) package from R. This version assumes all
    of the features are Gaussian.

    Read more in the :ref:`User Guide <mice>`.

    Parameters
    ----------
    missing_values : int, np.nan, optional (default=np.nan)
        The placeholder for the missing values. All occurrences of
        ``missing_values`` will be imputed.

    imputation_order : str, optional (default="ascending")
        The order in which the features will be imputed. Possible values:

        "ascending"
            From features with fewest missing values to most.
        "descending"
            From features with most missing values to fewest.
        "roman"
            Left to right.
        "arabic"
            Right to left.
        "random"
            A random order for each round.

    n_imputations : int, optional (default=100)
        Number of chained imputation rounds to perform, the results of which
        will be used in the final average.

    n_burn_in : int, optional (default=10)
        Number of initial imputation rounds to perform the results of which
        will not be returned.

    predictor : estimator object, default=BayesianRidge()
        The predictor to use at each step of the round-robin imputation.
        It must support ``return_std`` in its ``predict`` method.

    n_nearest_features : int, optional (default=None)
        Number of other features to use to estimate the missing values of
        the each feature column. Nearness between features is measured using
        the absolute correlation coefficient between each feature pair (after
        initial imputation). Can provide significant speed-up when the number
        of features is huge. If ``None``, all features will be used.

    initial_strategy : str, optional (default="mean")
        Which strategy to use to initialize the missing values. Same as the
        ``strategy`` parameter in :class:`sklearn.impute.SimpleImputer`
        Valid values: {"mean", "median", "most_frequent", or "constant"}.

    min_value : float, optional (default=None)
        Minimum possible imputed value. Default of ``None`` will set minimum
        to negative infinity.

    max_value : float, optional (default=None)
        Maximum possible imputed value. Default of ``None`` will set maximum
        to positive infinity.

    verbose : int, optional (default=0)
        Verbosity flag, controls the debug messages that are issued
        as functions are evaluated. The higher, the more verbose. Can be 0, 1,
        or 2.

    random_state : int, RandomState instance or None, optional (default=None)
        The seed of the pseudo random number generator to use when shuffling
        the data.  If int, random_state is the seed used by the random number
        generator; If RandomState instance, random_state is the random number
        generator; If None, the random number generator is the RandomState
        instance used by ``np.random``.

    Attributes
    ----------
    initial_imputer_ : object of class :class:`sklearn.preprocessing.Imputer`'
        The imputer used to initialize the missing values.

    imputation_sequence_ : list of tuples
        Each tuple has ``(feat_idx, neighbor_feat_idx, predictor)``, where
        ``feat_idx`` is the current feature to be imputed,
        ``neighbor_feat_idx`` is the array of other features used to impute the
        current feature, and ``predictor`` is the trained predictor used for
        the imputation.

    Notes
    -----
    The R version of MICE does not have inductive functionality, i.e. first
    fitting on ``X_train`` and then transforming any ``X_test`` without
    additional fitting. We do this by storing each feature's predictor during
    the round-robin ``fit`` phase, and predicting without refitting (in order)
    during the ``transform`` phase.

    Features which contain all missing values at ``fit`` are discarded upon
    ``transform``.

    Features with missing values in transform which did not have any missing
    values in fit will be imputed with the initial imputation method only.

    References
    ----------
    .. [1] `Stef van Buuren, Karin Groothuis-Oudshoorn (2011). "mice:
        Multivariate Imputation by Chained Equations in R". Journal of
        Statistical Software 45: 1-67.
        <https://www.jstatsoft.org/article/view/v045i03>`_
    """

    def __init__(self,
                 missing_values=np.nan,
                 imputation_order='ascending',
                 n_imputations=100,
                 n_burn_in=10,
                 predictor=None,
                 n_nearest_features=None,
                 initial_strategy="mean",
                 min_value=None,
                 max_value=None,
                 verbose=False,
                 random_state=None):

        self.missing_values = missing_values
        self.imputation_order = imputation_order
        self.n_imputations = n_imputations
        self.n_burn_in = n_burn_in
        self.predictor = predictor
        self.n_nearest_features = n_nearest_features
        self.initial_strategy = initial_strategy
        self.min_value = min_value
        self.max_value = max_value
        self.verbose = verbose
        self.random_state = random_state

    def _impute_one_feature(self,
                            X_filled,
                            mask_missing_values,
                            feat_idx,
                            neighbor_feat_idx,
                            predictor=None,
                            fit_mode=True):
        """Impute a single feature from the others provided.

        This function predicts the missing values of one of the features using
        the current estimates of all the other features. The ``predictor`` must
        support ``return_std=True`` in its ``predict`` method for this function
        to work.

        Parameters
        ----------
        X_filled : ndarray
            Input data with the most recent imputations.

        mask_missing_values : ndarray
            Input data's missing indicator matrix.

        feat_idx : int
            Index of the feature currently being imputed.

        neighbor_feat_idx : ndarray
            Indices of the features to be used in imputing ``feat_idx``.

        predictor : object
            The predictor to use at this step of the round-robin imputation.
            It must support ``return_std`` in its ``predict`` method.
            If None, it will be cloned from self._predictor.

        fit_mode : boolean, default=True
            Whether to fit and predict with the predictor or just predict.

        Returns
        -------
        X_filled : ndarray
            Input data with ``X_filled[missing_row_mask, feat_idx]`` updated.

        predictor : predictor with sklearn API
            The fitted predictor used to impute
            ``X_filled[missing_row_mask, feat_idx]``.
        """

        # if nothing is missing, just return the default
        # (should not happen at fit time because feat_ids would be excluded)
        missing_row_mask = mask_missing_values[:, feat_idx]
        if not np.any(missing_row_mask):
            return X_filled, predictor

        if predictor is None and fit_mode is False:
            raise ValueError("If fit_mode is False, then an already-fitted "
                             "predictor should be passed in.")

        if predictor is None:
            predictor = clone(self._predictor)

        if fit_mode:
            X_train = safe_indexing(X_filled[:, neighbor_feat_idx],
                                    ~missing_row_mask)
            y_train = safe_indexing(X_filled[:, feat_idx],
                                    ~missing_row_mask)
            predictor.fit(X_train, y_train)

        # get posterior samples
        X_test = safe_indexing(X_filled[:, neighbor_feat_idx],
                               missing_row_mask)
        mus, sigmas = predictor.predict(X_test, return_std=True)
        good_sigmas = sigmas > 0
        imputed_values = np.zeros(mus.shape, dtype=X_filled.dtype)
        imputed_values[~good_sigmas] = mus[~good_sigmas]
        imputed_values[good_sigmas] = self.random_state_.normal(
            loc=mus[good_sigmas], scale=sigmas[good_sigmas])

        # clip the values
        imputed_values = np.clip(imputed_values,
                                 self._min_value,
                                 self._max_value)

        # update the feature
        X_filled[missing_row_mask, feat_idx] = imputed_values
        return X_filled, predictor

    def _get_neighbor_feat_idx(self,
                               n_features,
                               feat_idx,
                               abs_corr_mat):
        """Get a list of other features to predict ``feat_idx``.

        If self.n_nearest_features is less than or equal to the total
        number of features, then use a probability proportional to the absolute
        correlation between ``feat_idx`` and each other feature to randomly
        choose a subsample of the other features (without replacement).

        Parameters
        ----------
        n_features : int
            Number of features in ``X``.

        feat_idx : int
            Index of the feature currently being imputed.

        abs_corr_mat : ndarray, shape (n_features, n_features)
            Absolute correlation matrix of ``X``. The diagonal has been zeroed
            out and each feature has been normalized to sum to 1. Can be None.

        Returns
        -------
        neighbor_feat_idx : array-like
            The features to use to impute ``feat_idx``.
        """
        if (self.n_nearest_features is not None and
                self.n_nearest_features < n_features):
            p = abs_corr_mat[:, feat_idx]
            neighbor_feat_idx = self.random_state_.choice(
                np.arange(n_features), self.n_nearest_features, replace=False,
                p=p)
        else:
            inds_left = np.arange(feat_idx)
            inds_right = np.arange(feat_idx + 1, n_features)
            neighbor_feat_idx = np.concatenate((inds_left, inds_right))
        return neighbor_feat_idx

    def _get_ordered_idx(self, mask_missing_values):
        """Decide in what order we will update the features.

        As a homage to the MICE R package, we will have 4 main options of
        how to order the updates, and use a random order if anything else
        is specified.

        Also, this function skips features which have no missing values.

        Parameters
        ----------
        mask_missing_values : array-like, shape (n_samples, n_features)
            Input data's missing indicator matrix, where "n_samples" is the
            number of samples and "n_features" is the number of features.

        Returns
        -------
        ordered_idx : ndarray, shape (n_features,)
            The order in which to impute the features.
        """
        frac_of_missing_values = mask_missing_values.mean(axis=0)
        missing_values_idx = np.nonzero(frac_of_missing_values)[0]
        if self.imputation_order == 'roman':
            ordered_idx = missing_values_idx
        elif self.imputation_order == 'arabic':
            ordered_idx = missing_values_idx[::-1]
        elif self.imputation_order == 'ascending':
            n = len(frac_of_missing_values) - len(missing_values_idx)
            ordered_idx = np.argsort(frac_of_missing_values,
                                     kind='mergesort')[n:][::-1]
        elif self.imputation_order == 'descending':
            n = len(frac_of_missing_values) - len(missing_values_idx)
            ordered_idx = np.argsort(frac_of_missing_values,
                                     kind='mergesort')[n:]
        elif self.imputation_order == 'random':
            ordered_idx = missing_values_idx
            self.random_state_.shuffle(ordered_idx)
        else:
            raise ValueError("Got an invalid imputation order: '{0}'. It must "
                             "be one of the following: 'roman', 'arabic', "
                             "'ascending', 'descending', or "
                             "'random'.".format(self.imputation_order))
        return ordered_idx

    def _get_abs_corr_mat(self, X_filled, tolerance=1e-6):
        """Get absolute correlation matrix between features.

        Parameters
        ----------
        X_filled : ndarray, shape (n_samples, n_features)
            Input data with the most recent imputations.

        tolerance : float, optional (default=1e-6)
            ``abs_corr_mat`` can have nans, which will be replaced
            with ``tolerance``.

        Returns
        -------
        abs_corr_mat : ndarray, shape (n_features, n_features)
            Absolute correlation matrix of ``X`` at the beginning of the
            current round. The diagonal has been zeroed out and each feature's
            absolute correlations with all others have been normalized to sum
            to 1.
        """
        n_features = X_filled.shape[1]
        if (self.n_nearest_features is None or
                self.n_nearest_features >= n_features):
            return None
        abs_corr_mat = np.abs(np.corrcoef(X_filled.T))
        # np.corrcoef is not defined for features with zero std
        abs_corr_mat[np.isnan(abs_corr_mat)] = tolerance
        # ensures exploration, i.e. at least some probability of sampling
        np.clip(abs_corr_mat, tolerance, None, out=abs_corr_mat)
        # features are not their own neighbors
        np.fill_diagonal(abs_corr_mat, 0)
        # needs to sum to 1 for np.random.choice sampling
        abs_corr_mat = normalize(abs_corr_mat, norm='l1', axis=0, copy=False)
        return abs_corr_mat

    def _initial_imputation(self, X):
        """Perform initial imputation for input X.

        Parameters
        ----------
        X : ndarray, shape (n_samples, n_features)
            Input data, where "n_samples" is the number of samples and
            "n_features" is the number of features.

        Returns
        -------
        Xt : ndarray, shape (n_samples, n_features)
            Input data, where "n_samples" is the number of samples and
            "n_features" is the number of features.

        X_filled : ndarray, shape (n_samples, n_features)
            Input data with the most recent imputations.

        mask_missing_values : ndarray, shape (n_samples, n_features)
            Input data's missing indicator matrix, where "n_samples" is the
            number of samples and "n_features" is the number of features.
        """
        if is_scalar_nan(self.missing_values):
            force_all_finite = "allow-nan"
        else:
            force_all_finite = True

        X = check_array(X, dtype=FLOAT_DTYPES, order="F",
                        force_all_finite=force_all_finite)
        _check_inputs_dtype(X, self.missing_values)

        mask_missing_values = _get_mask(X, self.missing_values)
        if self.initial_imputer_ is None:
            self.initial_imputer_ = SimpleImputer(
                                            missing_values=self.missing_values,
                                            strategy=self.initial_strategy)
            X_filled = self.initial_imputer_.fit_transform(X)
        else:
            X_filled = self.initial_imputer_.transform(X)

        valid_mask = np.flatnonzero(np.logical_not(
            np.isnan(self.initial_imputer_.statistics_)))
        Xt = X[:, valid_mask]
        mask_missing_values = mask_missing_values[:, valid_mask]

        return Xt, X_filled, mask_missing_values

    def fit_transform(self, X, y=None):
        """Fits the imputer on X and return the transformed X.

        Parameters
        ----------
        X : array-like, shape (n_samples, n_features)
            Input data, where "n_samples" is the number of samples and
            "n_features" is the number of features.

        y : ignored.

        Returns
        -------
        Xt : array-like, shape (n_samples, n_features)
             The imputed input data.
        """
        self.random_state_ = getattr(self, "random_state_",
                                     check_random_state(self.random_state))

        if self.predictor is None:
            from .linear_model import BayesianRidge
            self._predictor = BayesianRidge()
        else:
            self._predictor = clone(self.predictor)

        self._min_value = np.nan if self.min_value is None else self.min_value
        self._max_value = np.nan if self.max_value is None else self.max_value

        self.initial_imputer_ = None
        X, X_filled, mask_missing_values = self._initial_imputation(X)

        # edge case: in case the user specifies 0 for n_imputations,
        # then there is no need to do burn in and the result should be
        # just the initial imputation (before clipping)
        if self.n_imputations < 1:
            return X_filled

        X_filled = np.clip(X_filled, self._min_value, self._max_value)

        # order in which to impute
        # note this is probably too slow for large feature data (d > 100000)
        # and a better way would be good.
        # see: https://goo.gl/KyCNwj and subsequent comments
        ordered_idx = self._get_ordered_idx(mask_missing_values)

        abs_corr_mat = self._get_abs_corr_mat(X_filled)

        # impute data
        n_rounds = self.n_burn_in + self.n_imputations
        n_samples, n_features = X_filled.shape
        Xt = np.zeros((n_samples, n_features), dtype=X.dtype)
        self.imputation_sequence_ = []
        if self.verbose > 0:
            print("[ChainedImputer] Completing matrix with shape %s"
                  % (X.shape,))
        start_t = time()
        for i_rnd in range(n_rounds):
            if self.imputation_order == 'random':
                ordered_idx = self._get_ordered_idx(mask_missing_values)

            for feat_idx in ordered_idx:
                neighbor_feat_idx = self._get_neighbor_feat_idx(n_features,
                                                                feat_idx,
                                                                abs_corr_mat)
                X_filled, predictor = self._impute_one_feature(
                    X_filled, mask_missing_values, feat_idx, neighbor_feat_idx,
                    predictor=None, fit_mode=True)
                predictor_triplet = ImputerTriplet(feat_idx,
                                                   neighbor_feat_idx,
                                                   predictor)
                self.imputation_sequence_.append(predictor_triplet)

            if i_rnd >= self.n_burn_in:
                Xt += X_filled
            if self.verbose > 0:
                print('[ChainedImputer] Ending imputation round '
                      '%d/%d, elapsed time %0.2f'
                      % (i_rnd + 1, n_rounds, time() - start_t))

        Xt /= self.n_imputations
        Xt[~mask_missing_values] = X[~mask_missing_values]
        return Xt

    def transform(self, X):
        """Imputes all missing values in X.

        Note that this is stochastic, and that if random_state is not fixed,
        repeated calls, or permuted input, will yield different results.

        Parameters
        ----------
        X : array-like, shape = [n_samples, n_features]
            The input data to complete.

        Returns
        -------
        Xt : array-like, shape (n_samples, n_features)
             The imputed input data.
        """
        check_is_fitted(self, 'initial_imputer_')

        X, X_filled, mask_missing_values = self._initial_imputation(X)

        # edge case: in case the user specifies 0 for n_imputations,
        # then there is no need to do burn in and the result should be
        # just the initial imputation (before clipping)
        if self.n_imputations < 1:
            return X_filled

        X_filled = np.clip(X_filled, self._min_value, self._max_value)

        n_rounds = self.n_burn_in + self.n_imputations
        n_imputations = len(self.imputation_sequence_)
        imputations_per_round = n_imputations // n_rounds
        i_rnd = 0
        Xt = np.zeros(X.shape, dtype=X.dtype)
        if self.verbose > 0:
            print("[ChainedImputer] Completing matrix with shape %s"
                  % (X.shape,))
        start_t = time()
        for it, predictor_triplet in enumerate(self.imputation_sequence_):
            X_filled, _ = self._impute_one_feature(
                X_filled,
                mask_missing_values,
                predictor_triplet.feat_idx,
                predictor_triplet.neighbor_feat_idx,
                predictor=predictor_triplet.predictor,
                fit_mode=False
            )
            if not (it + 1) % imputations_per_round:
                if i_rnd >= self.n_burn_in:
                    Xt += X_filled
                if self.verbose > 1:
                    print('[ChainedImputer] Ending imputation round '
                          '%d/%d, elapsed time %0.2f'
                          % (i_rnd + 1, n_rounds, time() - start_t))
                i_rnd += 1

        Xt /= self.n_imputations
        Xt[~mask_missing_values] = X[~mask_missing_values]
        return Xt

    def fit(self, X, y=None):
        """Fits the imputer on X and return self.

        Parameters
        ----------
        X : array-like, shape (n_samples, n_features)
            Input data, where "n_samples" is the number of samples and
            "n_features" is the number of features.

        y : ignored

        Returns
        -------
        self : object
            Returns self.
        """
        self.fit_transform(X)
        return self


class MissingIndicator(BaseEstimator, TransformerMixin):
    """Binary indicators for missing values.

    Parameters
    ----------
    missing_values : number, string, np.nan (default) or None
        The placeholder for the missing values. All occurrences of
        `missing_values` will be imputed.

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
        in fit This is applicable only when ``features="missing-only"``.

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
    >>> indicator.fit(X1)
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

        if not ((isinstance(self.sparse, six.string_types) and
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
