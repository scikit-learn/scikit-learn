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
from ..neighbors.base import _get_weights, _check_weights
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
    """Imputation transformer for completing missing values using Nearest
    Neighbors. Broadly speaking, the imputation is performed using either
    the weighted or the unweighted mean of the desired number of neighbors.

    Parameters
    ----------
    missing_values : integer or "NaN", optional (default="NaN")
        The placeholder for the missing values. All occurrences of
        `missing_values` will be imputed. For missing values encoded as np.nan,
        use the string value "NaN".

    n_neighbors : int, optional (default = 5)
        Number of neighbors to get.

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
    >>> from sklearn.preprocessing import knn_imputation
    >>> nan = float("NaN")
    >>> X = [[1, 2, nan], [3, 4, 3], [nan, 6, 5], [8, 8, 7]]
    >>> imputer = knn_imputation.KNNImputer(n_neighbors=2, weights="uniform")
    >>> imputer.fit(X)
    KNNImputer(col_max_missing=0.8, copy=True, metric='masked_euclidean',
          missing_values='NaN', n_neighbors=2, row_max_missing=0.5,
          weights='uniform')
    >>> imputer.transform(X)
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
        self.weights = _check_weights(weights)
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
        # Check parameters
        force_all_finite = False if self.missing_values in ["NaN",
                                                            np.nan] else True
        X = check_array(X, accept_sparse=False, dtype=np.float64,
                        force_all_finite=force_all_finite, copy=self.copy)
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
        self.statistics_ = (neigh.fit(X), X_col_means)

        return self

    def transform(self, X):
        """Impute all missing values in X.

        Parameters
        ----------
        X : {array-like}, shape = [n_samples, n_features]
            The input data to complete.
        """
        force_all_finite = False if self.missing_values in ["NaN",
                                                            np.nan] else True
        check_is_fitted(self, 'statistics_')
        X = check_array(X, accept_sparse=False, dtype=FLOAT_DTYPES,
                        force_all_finite=force_all_finite, copy=self.copy)
        mask = _get_mask(X, self.missing_values)
        n_rows_X, n_cols_X = X.shape
        row_total_missing = mask.sum(axis=1)

        # Get fitted objects
        fitted_NNObj, fitted_col_means = self.statistics_
        fitted_X = fitted_NNObj._fit_X

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
        if X is fitted_X or X.base is fitted_X.base:
            neighbors = fitted_NNObj.kneighbors(n_neighbors=self.n_neighbors)
        else:
            neighbors = fitted_NNObj.kneighbors(X,
                                                n_neighbors=self.n_neighbors)

        # Get row index, distance, and weights of donors
        knn_distances, knn_row_index = neighbors
        weights = _get_weights(knn_distances[row_total_missing.astype(
            np.bool), ], self.weights)

        knn_row_index = np.vsplit(knn_row_index, knn_row_index.shape[0])
        knn_row_index = np.repeat(knn_row_index,
                                  row_total_missing, axis=0).ravel()

        # Get column index of donors
        # NOTE: Following assumes columns in X and _fit_X are in the same order
        col_missing_index = np.where(mask)[1]
        knn_col_index = np.repeat(col_missing_index, self.n_neighbors)

        # Calculate kNN score and impute
        donors = fitted_X[
            (knn_row_index, knn_col_index)].reshape((-1, self.n_neighbors))
        donors = np.ma.array(
            donors, mask=_get_mask(donors, self.missing_values))
        imputed = np.ma.average(donors, axis=1, weights=weights)
        X[mask] = imputed.data

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
