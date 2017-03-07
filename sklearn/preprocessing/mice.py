# Authors: Sergey Feldman <sergeyfeldman@gmail.com>
# License: BSD 3 clause

from time import time

import numpy as np

from .imputation import _get_mask, Imputer
from ..base import BaseEstimator, TransformerMixin, clone
from ..dummy import DummyRegressor
from ..preprocessing import normalize
from ..utils import check_array, check_random_state
from ..utils.random import choice
from ..utils.validation import check_is_fitted

__all__ = [
    'MICEImputer',
]


class MICEImputer(BaseEstimator, TransformerMixin):
    """Basic implementation of MICE package from R.
    This version assumes all of the columns are Gaussian.

    Parameters
    ----------
    missing_values : int or "NaN", optional (default="NaN")
        The placeholder for the missing values. All occurrences of
        missing_values will be imputed. For missing values encoded as
        np.nan, use the string value "NaN".

    imputation_order : str, optional (default="monotone")
        The order in which the columns will be imputed.

        - "monotone" - From columns with fewest missing values to most.
        - "revmonotone" - From columns with most missing values to fewest.
        - "roman" - Left to right.
        - "arabic" - Right to left.
        - "random" - A random order for each round.

    n_imputations : int, optional (defaults = 100)
        Number of MICE rounds to perform the results of which will be
        used in the final average.

    n_burn_in : int, optional (default = 10)
        Number of initial MICE rounds to perform the results of which
        will not be returned.

    n_nearest_columns : int, optional (default = np.infty)
        Number of other columns to use to estimate the missing values of
        the current column. Can provide significant speed-up
        when the number of columns is huge.

    initial_fill_method : str, optional (default = "mean" )
        Valid values: {"mean", "median", or "most_frequent"}.
        Uses :class:sklearn.preprocessing.Imputer.

    min_value : float (default = -np.infty)
        Minimum possible imputed value.

    max_value : float (default = np.infty)
        Maximum possible imputed value.

    verbose : boolean (default = False)
        Whether to provide updates on progress.

    random_state : int seed, RandomState instance (default = None)
        The seed of the pseudo random number generator to use when
        shuffling the data.

    Notes
    -----
    Columns which only contain missing values at `fit` are discarded upon
        `transform`.

    References
    ----------
    .. [1] `Stef van Buuren, Karin Groothuis-Oudshoorn (2011). "mice:
        Multivariate Imputation by Chained Equations in R". Journal of
        Statistical Software 45: 1-67.
        <https://www.jstatsoft.org/article/view/v045i03>`_
    """

    def __init__(
            self,
            missing_values='NaN',
            imputation_order='monotone',
            n_imputations=100,
            n_burn_in=10,
            n_nearest_columns=np.infty,
            initial_fill_method="mean",
            min_value=-np.infty,
            max_value=np.infty,
            verbose=False,
            random_state=None):
        from ..linear_model import BayesianRidge  # avoiding circular import
        self.model = BayesianRidge()
        self.missing_values = missing_values
        self.imputation_order = imputation_order
        self.n_imputations = n_imputations
        self.n_burn_in = n_burn_in
        self.n_nearest_columns = n_nearest_columns
        self.initial_fill_method = initial_fill_method
        self.min_value = min_value
        self.max_value = max_value
        self.verbose = verbose
        self.random_state = random_state

    def _fill_in_one_column(self,
                            X_filled,
                            mask_missing_values,
                            this_column,
                            other_columns,
                            model=None,
                            min_std=1e-6):
        """
        This function predicts the missing values of one of the
        features using the current estimates of all the other
        features. The `model` must support `return_std=True`
        in its `predict` method for this function to work.
        """

        # if nothing is missing, just return the default
        rng = check_random_state(self.random_state)
        if mask_missing_values[:, this_column].sum() == 0:
            return X_filled, model
        missing_row_mask = mask_missing_values[:, this_column]
        if model is None:
            X_train = X_filled[:, other_columns][~missing_row_mask]
            y_train = X_filled[:, this_column][~missing_row_mask]
            if np.std(y_train) > 0:
                model = clone(self.model)
                model.fit(X_train, y_train)
            else:
                model = DummyRegressor()
                model.fit(X_train, y_train)
        X_test = X_filled[:, other_columns][missing_row_mask]
        mus, sigmas = model.predict(X_test, return_std=True)
        sigmas = np.maximum(sigmas, min_std)
        imputed_values = rng.normal(mus, sigmas)
        imputed_values = np.clip(imputed_values,
                                 self.min_value,
                                 self.max_value)
        X_filled[missing_row_mask, this_column] = imputed_values
        return X_filled, model

    def _get_other_columns(self,
                           n_features,
                           this_column,
                           abs_correlation_matrix):
        """Gets a list of other columns to predict this_column.

        If self.n_nearest_columns is less than or equal to the total
        number of columns, then use a probability proportional to the absolute
        correlation between this_column and each other column to randomly
        choose a subsample of the other columns (without replacement).

        Parameters
        ----------
        n_features : integer
            Number of columns in X.

        this_column : integer
            Index of the column currently being imputed.

        abs_correlation_matrix : array-like, shape (n_features, n_features)
            Absolute correlation matrix of X at the begging of the current
            round. The diagonal has been zeroed out and each column has been
            normalized to sum to 1.
        """
        ordered_columns = np.arange(n_features)
        if self.n_nearest_columns <= n_features - 1:
            p = abs_correlation_matrix[:, this_column]
            other_columns = choice(ordered_columns,
                                   self.n_nearest_columns,
                                   replace=False,
                                   p=p)
        else:
            other_columns = np.concatenate([ordered_columns[:this_column],
                                            ordered_columns[this_column + 1:]])
        return other_columns

    def _get_ordered_indices(self, mask_missing_values):
        """Decides in what order we will update the columns.

        As a homage to the MICE R package, we will have 4 main options of
        how to order the updates, and use a random order if anything else
        is specified.

        Also, this function filters out columns which have no missing values.

        Parameters
        ----------
        mask_missing_values : array-like, shape (n_samples, n_features)
            Input data's missing indicator matrix, where "n_samples" is the
            number of samples and "n_features" is the number of features.
        """
        rng = check_random_state(self.random_state)
        n_samples, n_features = mask_missing_values.shape
        fraction_of_missing_values = mask_missing_values.mean(axis=0)
        all_features_indices = np.arange(n_features)
        if self.imputation_order == 'roman':
            ordered_indices = all_features_indices
        elif self.imputation_order == 'arabic':
            ordered_indices = all_features_indices[::-1]
        elif self.imputation_order == 'monotone':
            ordered_indices = np.argsort(fraction_of_missing_values)[::-1]
        elif self.imputation_order == 'revmonotone':
            ordered_indices = np.argsort(fraction_of_missing_values)
        else:
            ordered_indices = np.arange(n_features)
            rng.shuffle(ordered_indices)

        # filter out indices for which we have no missing values
        valid_features = all_features_indices[fraction_of_missing_values > 0]
        valid_features = set(valid_features)
        ordered_indices = [i for i in ordered_indices if i in valid_features]
        return ordered_indices

    def _get_abs_correlation_matrix(self, X_filled, eps=1e-6):
        # at each stage all but one of the features is used as input
        n_features = X_filled.shape[1]
        if self.n_nearest_columns > n_features - 1:
            return None
        abs_correlation_matrix = np.abs(np.corrcoef(X_filled.T))
        # np.corrcoef is not defined for constant columns
        abs_correlation_matrix[np.isnan(abs_correlation_matrix)] = eps
        np.fill_diagonal(abs_correlation_matrix, 0)
        abs_correlation_matrix = normalize(abs_correlation_matrix, norm='l1',
                                           axis=0)
        return abs_correlation_matrix

    def fit_transform(self, X, y=None):
        """Fits the imputer on X and return the transformed X.

        Parameters
        ----------
        X : array-like, shape (n_samples, n_features)
            Input data, where "n_samples" is the number of samples and
            "n_features" is the number of features.
        """
        X = check_array(X, dtype=np.float64, force_all_finite=False)
        X = np.asarray(X, order="F")
        mask_missing_values = _get_mask(X, self.missing_values)
        self.trained_model_triplets = []

        # initial imputation
        self.initial_imputer_ = Imputer(missing_values=self.missing_values,
                                        strategy=self.initial_fill_method,
                                        axis=0)
        X_filled = self.initial_imputer_.fit_transform(X)
        self._val_inds = self.initial_imputer_._valid_statistics_indexes
        X = X[:, self._val_inds]
        mask_missing_values = mask_missing_values[:, self._val_inds]

        # perform imputations
        n_samples, n_features = X_filled.shape
        total_rounds = self.n_burn_in + self.n_imputations
        results_list = []
        if self.verbose:
            print("[MICE] Completing matrix with shape %s" % (X.shape,))
            start_t = time()
            mice_msg = '[MICE] Ending imputation round '
        for m in range(total_rounds):
            # order in which to impute
            ordered_indices = self._get_ordered_indices(mask_missing_values)

            # abs_correlation matrix is used to choose a subset of other
            # features to impute from  f
            abs_corr_mat = self._get_abs_correlation_matrix(X_filled)

            # Fill in each column in the order of ordered_indices
            for this_column in ordered_indices:
                other_columns = self._get_other_columns(n_features,
                                                        this_column,
                                                        abs_corr_mat)
                X_filled, model = self._fill_in_one_column(X_filled,
                                                           mask_missing_values,
                                                           this_column,
                                                           other_columns)
                model_triplet = (this_column, other_columns, model)
                self.trained_model_triplets.append(model_triplet)

            if m >= self.n_burn_in:
                results_list.append(X_filled[mask_missing_values])
            if self.verbose:
                print(mice_msg + 'round %d/%d, elapsed time %0.2f'
                      % (m + 1, total_rounds, time() - start_t))

        if len(results_list) > 0:
            X[mask_missing_values] = np.array(results_list).mean(axis=0)
        else:
            X[mask_missing_values] = X_filled[mask_missing_values]

        return X

    def transform(self, X):
        """Imputes all missing values in X.

        Parameters
        ----------
        X : array-like}, shape = [n_samples, n_features]
            The input data to complete.
        """
        check_is_fitted(self, 'initial_imputer_')
        X = check_array(X, dtype=np.float64, force_all_finite=False)
        X = np.asarray(X, order="F")
        mask_missing_values = _get_mask(X, self.missing_values)

        # initial imputation
        X_filled = self.initial_imputer_.transform(X)
        X = X[:, self._val_inds]
        mask_missing_values = mask_missing_values[:, self._val_inds]

        # perform imputations
        n_samples, n_features = X_filled.shape
        total_rounds = self.n_burn_in + self.n_imputations
        if total_rounds > 0:
            total_iterations = len(self.trained_model_triplets)
            imputations_per_round = total_iterations / total_rounds
            round_index = 0
            results_list = []
            if self.verbose:
                print("[MICE] Completing matrix with shape %s" % (X.shape,))
                start_t = time()
                mice_msg = '[MICE] Ending imputation round '
            for i, model_triplet in enumerate(self.trained_model_triplets):
                this_column, other_columns, model = model_triplet
                X_filled, _ = self._fill_in_one_column(X_filled,
                                                       mask_missing_values,
                                                       this_column,
                                                       other_columns,
                                                       model)
                if not (i + 1) % imputations_per_round:
                    round_index += 1
                    if round_index >= self.n_burn_in:
                        results_list.append(X_filled[mask_missing_values])
                    if self.verbose:
                        print(mice_msg + '%d/%d, elapsed time %0.2f'
                              % (round_index, total_rounds, time() - start_t))

        if total_rounds > 0 and len(results_list) > 0:
            X[mask_missing_values] = np.array(results_list).mean(axis=0)
        else:
            X[mask_missing_values] = X_filled[mask_missing_values]

        return X

    def fit(self, X, y=None):
        """Fits the imputer on X and return self.

        Parameters
        ----------
        X : array-like, shape (n_samples, n_features)
            Input data, where "n_samples" is the number of samples and
            "n_features" is the number of features.
        """
        self.fit_transform(X)
        return self
