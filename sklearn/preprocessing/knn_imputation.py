# Authors: Ashim Bhattarai <ashimb9@gmail.com>
# License: BSD 3 clause

from __future__ import division
import warnings
import numpy as np

from ..base import BaseEstimator, TransformerMixin
from ..utils import check_array
from ..utils.validation import check_is_fitted
from ..utils.validation import FLOAT_DTYPES
from ..neighbors import NearestNeighbors

__all__ = [
    'KNNImputer',
]


def _get_mask(X, value_to_mask):
    """Compute the boolean mask X == missing_values."""
    if value_to_mask == "NaN" or np.isnan(value_to_mask):
        return np.isnan(X)
    else:
        return X == value_to_mask


class KNNImputer(BaseEstimator, TransformerMixin):
    """Imputation transformer for completing missing values.

    Read more in the :ref:`User Guide <imputation>`.

    Parameters
    ----------
    missing_values : integer or "NaN", optional (default="NaN")
        The placeholder for the missing values. All occurrences of
        `missing_values` will be imputed. For missing values encoded as np.nan,
        use the string value "NaN".

    n_neighbors : int, optional (default = 10)
        Number of neighbors to get.

    weighted : bool, optional (default = True)
        Set the imputed value as a distance-weighted score of the neighbors

    metric : string or callable, optional (default = 'masked_euclidean')
        metric to use for distance computation.

    row_max_missing: float, optional (default = 0.5)
        The maximum percentage of columns (i.e. features) that can be missing
        before the sample is excluded from nearest neighbor imputation. It
        means that such rows will not be considered a potential donor in fit()
        and in transform() their missing feature values will be imputed to be
        the column mean for the entire dataset.

    col_max_missing: float, optional (default = 0.8)
        The maximum percentage of rows (or samples) that can be missing
        for a given feature beyond which an error is raised.

    copy : boolean, optional (default=True)
        If True, a copy of X will be created. If False, imputation will
        be done in-place whenever possible. Note that, if metric is
        "masked_euclidean" and copy=False then missing_values in the
        input matrix X will be overwritten with zeros.

    Attributes
    ----------
    statistics_ : array of shape (n_features,)
        A tuple whose first element if the fitted NearestNeighbors object
        and second element is the column means using available values.

    Notes
    -----
    """

    def __init__(self, missing_values="NaN", n_neighbors=10,
                 weighted=True, metric="masked_euclidean",
                 row_max_missing=0.5, col_max_missing=0.8, copy=True):

        self.missing_values = missing_values
        self.n_neighbors = n_neighbors
        self.weighted = weighted
        self.metric = metric
        self.row_max_missing = row_max_missing
        self.col_max_missing = col_max_missing
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
        X = check_array(X, accept_sparse=False, dtype=np.float64,
                        force_all_finite=False, copy=self.copy)
        mask = _get_mask(X, self.missing_values)

        # Check if % missing in any column > col_max_missing
        if np.any(mask.sum(axis=0) > (X.shape[0] * self.col_max_missing)):
            raise ValueError("The following columns have, "
                             "more than {0}% missing values: {1}"
                             .format(self.col_max_missing*100, np.where(
                              mask.sum(axis=0) > (X.shape[0] *
                                                  self.col_max_missing))))
        # X_masked = np.ma.array(X, mask=mask)
        X_col_means = X.mean(axis=0)

        # Check if % missing in any row > col_max_missing
        bad_rows = mask.sum(axis=1) > (mask.shape[1] * self.row_max_missing)
        if np.any(bad_rows):
            warnings.warn(
                "The following rows have more than {0}% missing values and "
                "are not included as nearest neighbors: {1}"
                    .format(self.row_max_missing*100, np.where(bad_rows)))

            # Remove rows that have more than row_max_missing % missing
            X = X[~bad_rows, :]
            mask = _get_mask(X, self.missing_values)
            # X_masked = np.ma.array(X, mask=mask)

        if X.shape[0] < self.n_neighbors:
            raise ValueError("There are only %d samples, "
                             "but n_neighbors=%d."
                             % (X.shape[0], self.n_neighbors))

        # Instantiate NN object, get column means, and store in statistics_
        neigh = NearestNeighbors(n_neighbors=self.n_neighbors,
                                 metric=self.metric)
        self.statistics_ = (neigh.fit(X), X_col_means)

        return self

    def transform(self, X):
        """Impute all missing values in X.

        Parameters
        ----------
        X : {array-like, sparse matrix}, shape = [n_samples, n_features]
            The input data to complete.
        """
        check_is_fitted(self, 'statistics_')
        X = check_array(X, accept_sparse=False, dtype=FLOAT_DTYPES,
                        force_all_finite=False, copy=self.copy)
        mask = _get_mask(X, self.missing_values)
        n_rows_X, n_cols_X = X.shape
        row_total_missing = mask.sum(axis=1)

        # Get fitted objects
        fitted_NNObj, fitted_col_means = self.statistics_
        fitted_X = fitted_NNObj._fit_X
        fitted_mask = _get_mask(fitted_X, self.missing_values)

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

        # Check if the X in fit() and transform() are the same
        if np.ma.allequal(np.ma.array(X, mask=mask),
                          np.ma.array(fitted_X, mask=fitted_mask)) and \
                np.array_equal(mask, fitted_mask):

            # Get the k nearest neighbors from fitted matrix
            neighbors = fitted_NNObj.kneighbors(n_neighbors=self.n_neighbors,
                                                return_distance=self.weighted)
        else:
            neighbors = fitted_NNObj.kneighbors(X,
                                                n_neighbors=self.n_neighbors,
                                                return_distance=self.weighted)

        # Get row index and distance (if weighted) of donors
        if self.weighted:
            knn_row_index = neighbors[1]
            knn_distances = neighbors[0]
        else:
            knn_row_index = neighbors
            # knn_distances = np.ones_like(neighbors)

        knn_row_index = np.vsplit(knn_row_index, knn_row_index.shape[0])
        knn_row_index = np.repeat(knn_row_index, row_total_missing, axis=0)
        knn_row_index = knn_row_index.ravel()

        # Get column index of donors
        # NOTE: Following assumes columns in X and _fit_X are in the same order
        col_missing_index = np.where(mask)[1]
        knn_col_index = np.repeat(col_missing_index, self.n_neighbors)

        # Calculate kNN score and impute
        imputed = np.nanmean(
            (fitted_NNObj._fit_X[(knn_row_index, knn_col_index)]).
                             reshape((-1, self.n_neighbors)), axis=1)
        X[mask] = imputed

        # Merge bad rows to X and mean impute any leftover missing
        if np.any(bad_rows):
            X_merged = np.empty((n_rows_X, n_cols_X))
            X_merged[bad_rows, :] = X_bad
            X_merged[~bad_rows, :] = X
            X = X_merged

        # Impute bad_rows and leftover missing with column means
        mask_after_knn = _get_mask(X, self.missing_values)
        if np.any(mask_after_knn):
            missing_index = np.where(mask_after_knn)
            X[missing_index] = np.take(fitted_col_means, missing_index[1])

        return X