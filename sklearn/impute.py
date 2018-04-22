"""Transformers for missing value imputation"""
# Authors: Nicolas Tresegnie <nicolas.tresegnie@gmail.com>
#          Sergey Feldman <sergeyfeldman@gmail.com>
#          Ashim Bhattarai <"ashimb9" + "\100gmail\56com">
# License: BSD 3 clause

from __future__ import division

import warnings
from time import time

import numpy as np
import numpy.ma as ma
from scipy import sparse
from scipy import stats
from collections import namedtuple

from .base import BaseEstimator, TransformerMixin
from .utils import check_array, _get_n_jobs
from .base import clone
from .preprocessing import normalize
from .utils import check_array, check_random_state, safe_indexing
from .utils.sparsefuncs import _get_median
from .utils.validation import check_is_fitted
from .utils.validation import FLOAT_DTYPES
from .metrics import pairwise_distances

from .neighbors.base import _check_weights
from .neighbors.base import _get_weights

from .externals import six

zip = six.moves.zip
map = six.moves.map

MICETriplet = namedtuple('MICETriplet', ['feat_idx',
                                         'neighbor_feat_idx',
                                         'predictor'])

__all__ = [
    'SimpleImputer',
    'MICEImputer',
    'KNNImputer',
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


class SimpleImputer(BaseEstimator, TransformerMixin):
    """Imputation transformer for completing missing values.

    Read more in the :ref:`User Guide <impute>`.

    Parameters
    ----------
    missing_values : integer or "NaN", optional (default="NaN")
        The placeholder for the missing values. All occurrences of
        `missing_values` will be imputed. For missing values encoded as np.nan,
        use the string value "NaN".

    strategy : string, optional (default="mean")
        The imputation strategy.

        - If "mean", then replace missing values using the mean along
          each column.
        - If "median", then replace missing values using the median along
          each column.
        - If "most_frequent", then replace missing using the most frequent
          value along each column.

    verbose : integer, optional (default=0)
        Controls the verbosity of the imputer.

    copy : boolean, optional (default=True)
        If True, a copy of X will be created. If False, imputation will
        be done in-place whenever possible. Note that, in the following cases,
        a new copy will always be made, even if `copy=False`:

        - If X is not an array of floating values;
        - If X is sparse and `missing_values=0`;
        - If X is encoded as a CSR matrix.

    Attributes
    ----------
    statistics_ : array of shape (n_features,)
        The imputation fill value for each feature.

    Notes
    -----
    Columns which only contained missing values at `fit` are discarded upon
    `transform`.

    """
    def __init__(self, missing_values="NaN", strategy="mean",
                 verbose=0, copy=True):
        self.missing_values = missing_values
        self.strategy = strategy
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
        self : SimpleImputer
        """
        # Check parameters
        allowed_strategies = ["mean", "median", "most_frequent"]
        if self.strategy not in allowed_strategies:
            raise ValueError("Can only use these strategies: {0} "
                             " got strategy={1}".format(allowed_strategies,
                                                        self.strategy))

        X = check_array(X, accept_sparse='csc', dtype=FLOAT_DTYPES,
                        force_all_finite='allow-nan'
                        if self.missing_values == 'NaN'
                        or np.isnan(self.missing_values) else True)

        if sparse.issparse(X):
            self.statistics_ = self._sparse_fit(X,
                                                self.strategy,
                                                self.missing_values)
        else:
            self.statistics_ = self._dense_fit(X,
                                               self.strategy,
                                               self.missing_values)

        return self

    def _sparse_fit(self, X, strategy, missing_values):
        """Fit the transformer on sparse data."""
        # Count the zeros
        if missing_values == 0:
            n_zeros_axis = np.zeros(X.shape[1], dtype=int)
        else:
            n_zeros_axis = X.shape[0] - np.diff(X.indptr)

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
                sums = X.sum(axis=0)
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

    def _dense_fit(self, X, strategy, missing_values):
        """Fit the transformer on dense data."""
        X = check_array(X, force_all_finite='allow-nan'
                        if self.missing_values == 'NaN'
                        or np.isnan(self.missing_values) else True)
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
        check_is_fitted(self, 'statistics_')
        X = check_array(X, accept_sparse='csc', dtype=FLOAT_DTYPES,
                        force_all_finite='allow-nan'
                        if self.missing_values == 'NaN'
                        or np.isnan(self.missing_values) else True,
                        copy=self.copy)
        statistics = self.statistics_
        if X.shape[1] != statistics.shape[0]:
            raise ValueError("X has %d features per sample, expected %d"
                             % (X.shape[1], self.statistics_.shape[0]))

        # Delete the invalid columns
        invalid_mask = np.isnan(statistics)
        valid_mask = np.logical_not(invalid_mask)
        valid_statistics = statistics[valid_mask]
        valid_statistics_indexes = np.flatnonzero(valid_mask)
        missing = np.arange(X.shape[1])[invalid_mask]

        if invalid_mask.any():
            if self.verbose:
                warnings.warn("Deleting features without "
                              "observed values: %s" % missing)
            X = X[:, valid_statistics_indexes]

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
            n_missing = np.sum(mask, axis=0)
            values = np.repeat(valid_statistics, n_missing)

            coordinates = np.where(mask.transpose())[::-1]

            X[coordinates] = values

        return X



class MICEImputer(BaseEstimator, TransformerMixin):
    """MICE transformer to impute missing values.

    Basic implementation of MICE (Multivariate Imputations by Chained
    Equations) package from R. This version assumes all of the features are
    Gaussian.

    Read more in the :ref:`User Guide <mice>`.

    Parameters
    ----------
    missing_values : int or "NaN", optional (default="NaN")
        The placeholder for the missing values. All occurrences of
        ``missing_values`` will be imputed. For missing values encoded as
        np.nan, use the string value "NaN".

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
        Number of MICE rounds to perform, the results of which will be
        used in the final average.

    n_burn_in : int, optional (default=10)
        Number of initial MICE rounds to perform the results of which
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
        ``strategy`` parameter in :class:`sklearn.preprocessing.Imputer`
        Valid values: {"mean", "median", or "most_frequent"}.

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
                 missing_values='NaN',
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
        abs_corr_mat[abs_corr_mat < tolerance] = tolerance
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
        X = check_array(X, dtype=FLOAT_DTYPES, order="F",
                        force_all_finite='allow-nan'
                        if self.missing_values == 'NaN'
                        or np.isnan(self.missing_values) else True)

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
            print("[MICE] Completing matrix with shape %s" % (X.shape,))
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
                predictor_triplet = MICETriplet(feat_idx,
                                                neighbor_feat_idx,
                                                predictor)
                self.imputation_sequence_.append(predictor_triplet)

            if i_rnd >= self.n_burn_in:
                Xt += X_filled
            if self.verbose > 0:
                print('[MICE] Ending imputation round '
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
            print("[MICE] Completing matrix with shape %s" % (X.shape,))
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
                    print('[MICE] Ending imputation round '
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
    missing_values : integer or "NaN", optional (default = "NaN")
        The placeholder for the missing values. All occurrences of
        `missing_values` will be imputed. For missing values encoded as
        ``np.nan``, use the string value "NaN".

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
            # Argpartition to seperate actual donors from the rest
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

            # Retrieve donor cells and calculate kNN score
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
        X : {array-like}, shape (n_samples, n_features)
            Input data, where ``n_samples`` is the number of samples and
            ``n_features`` is the number of features.

        Returns
        -------
        self : object
            Returns self.
        """

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
        X : {array-like}, shape = [n_samples, n_features]
            The imputed dataset.
        """

        check_is_fitted(self, ["fitted_X_", "statistics_"])
        force_all_finite = False if self.missing_values in ["NaN",
                                                            np.nan] else True
        X = check_array(X, accept_sparse=False, dtype=FLOAT_DTYPES,
                        force_all_finite=force_all_finite, copy=self.copy)
        # Check for +/- inf
        if np.any(np.isinf(X)):
            raise ValueError("+/- inf values are not allowed in data to be "
                             "transformed.")

        # Get fitted data and ensure correct dimension
        # fitted_X = self._fitted_neighbors._fit_X
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

            # Get row index of missing and distance from donors
            dist_temp = pairwise_distances(X[row_has_missing],
                                           self.fitted_X_,
                                           metric=self.metric,
                                           squared=False,
                                           missing_values=self.missing_values)
            dist = np.empty((n_rows_X, n_rows_X))
            dist[row_has_missing] = dist_temp.copy()
            dist[~row_has_missing] = np.nan

            # Delete temp var binding
            del dist_temp

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
        X : {array-like}, shape (n_samples, n_features)
            Input data, where ``n_samples`` is the number of samples and
            ``n_features`` is the number of features.

        Returns
        -------
        X : {array-like}, shape (n_samples, n_features)
            Returns imputed dataset.
        """
        return self.fit(X).transform(X)
