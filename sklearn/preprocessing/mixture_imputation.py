# Authors: Ashim Bhattarai <ashimb9@gmail.com>
# License: BSD 3 clause

import warnings

import numpy as np
# import numpy.ma as ma
from scipy import sparse
# from scipy import stats

from ..base import BaseEstimator, TransformerMixin
from ..utils import check_array
# from ..utils.sparsefuncs import _get_median
from ..utils.validation import check_is_fitted
from ..utils.validation import FLOAT_DTYPES
from ..mixture import gaussian_mixture

# from ..externals import six

# zip = six.moves.zip
# map = six.moves.map

__all__ = [
    'MixtureImputer',
]


def _get_mask(X, value_to_mask):
    """Compute the boolean mask X == missing_values."""
    if value_to_mask == "NaN" or np.isnan(value_to_mask):
        return np.isnan(X)
    else:
        return X == value_to_mask


class MixtureImputer(BaseEstimator, TransformerMixin):
    def __init__(self, missing_values="NaN", dist="gaussian",
                 axis=0, verbose=1, copy=True, min_complete=10,
                 n_components=1, covariance_type="full",
                 tol=0.001, max_iter=100):

        self.missing_values = missing_values
        self.dist = dist
        self.axis = axis
        self.verbose = verbose
        self.copy = copy
        self.min_complete = min_complete
        self.n_components = n_components
        self.covariance_type = covariance_type
        self.tol = tol
        self.max_iter = max_iter

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
                raise ValueError("Sparse matrices not supported yet")
                # self.statistics_ = self._sparse_fit(X,
                #                                     self.dist,
                #                                     self.missing_values,
                #                                     self.axis)
            else:
                self.statistics_ = self._dense_fit(X,
                                                   self.dist,
                                                   self.missing_values,
                                                   self.axis)

        return self

    def _dense_fit(self, X, dist, missing_values, axis):
        """Fit the transformer on dense data."""
        X = check_array(X, force_all_finite=False)
        if(self.copy):
            X = np.copy(X)
        mask = _get_mask(X, missing_values)
        # masked_X = ma.masked_array(X, mask=mask)
        other_axis = int(not axis)
        is_missing_sample = (mask.sum(axis=other_axis) > 0)
        n_complete_samples = np.sum(~is_missing_sample)

        if n_complete_samples < self.min_complete:
            raise ValueError("Must have {0} complete samples, "
                             " but total complete samples={1}".
                             format(self.min_complete, n_complete_samples))

        imputer = {"gaussian": gaussian_mixture.GaussianMixture
                   }[dist](n_components=self.n_components,
                                covariance_type=self.covariance_type,
                                tol=self.tol, max_iter=self.max_iter)

        # Initial fit with complete cases
        imputer.fit(X[~is_missing_sample, :])

        return imputer

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
            # if X.shape[1] != statistics.shape[0]:
            #     raise ValueError("X has %d features per sample, expected %d"
            #                      % (X.shape[1], self.statistics_.shape[0]))

        # Since two different arrays can be provided in fit(X) and
        # transform(X), the imputation data need to be recomputed
        # when the imputation is done per sample
        else:
            X = check_array(X, accept_sparse='csr', dtype=FLOAT_DTYPES,
                            force_all_finite=False, copy=self.copy)

            if sparse.issparse(X):
                raise ValueError("Sparse matrices not supported yet")

            else:
                statistics = self._dense_fit(X,
                                             self.strategy,
                                             self.missing_values,
                                             self.axis)

        mixture_means = statistics.means_
        # mixture_weights = statistics.weights_
        means_weighted_log_prob = statistics. \
            _estimate_weighted_log_prob(mixture_means)

        # The weighted mean value
        prediction = np.sum(np.exp(means_weighted_log_prob) *
                            mixture_means) / np.exp(statistics.
                                                    score_samples(mixture_means))

        # Delete the invalid rows/columns
        invalid_mask = np.isnan(prediction)
        valid_mask = np.logical_not(invalid_mask)
        valid_statistics = prediction[valid_mask]
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

        # Get replancement index and impute initial prediction
        mask = _get_mask(X, self.missing_values)
        replace_indices = np.where(mask)[not self.axis]
        X[mask] = np.take(valid_statistics, replace_indices)

        # Fit then impute (repeat until convergence)
        for _ in range(self.max_iter):
            loglik_old = statistics.score(X)
            statistics.fit(X)
            mixture_means = statistics.means_
            # mixture_weights = statistics.weights_
            means_weighted_log_prob = statistics. \
                _estimate_weighted_log_prob(mixture_means)

            # The weighted mean value
            prediction = np.sum(np.exp(means_weighted_log_prob) *
                                mixture_means) / np.exp(statistics.
                                                        score_samples(mixture_means))
            # Impute
            # replace_indices = np.where(mask)[not self.axis]
            X[mask] = np.take(prediction, replace_indices)

            loglik_new = statistics.score(X)
            if (loglik_new - loglik_old) < self.tol:
                break

        return X