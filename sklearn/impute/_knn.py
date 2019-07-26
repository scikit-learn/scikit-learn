import numpy as np

from ..base import BaseEstimator, TransformerMixin
from ..utils.validation import FLOAT_DTYPES
from ..metrics import pairwise_distances
from ..metrics.pairwise import _NAN_METRICS
from ..neighbors.base import _get_weights
from ..neighbors.base import _check_weights
from ..utils import check_array
from ..utils import is_scalar_nan
from ..utils.mask import _get_missing_mask
from ..utils.validation import check_is_fitted


class KNNImputer(TransformerMixin, BaseEstimator):
    """Imputation for completing missing values using k-Nearest Neighbors.

    Each sample's missing values are imputed using values from ``n_neighbors``
    nearest neighbors found in the training set. Each missing feature is then
    imputed as the average, either weighted or unweighted, of these neighbors.
    If a sample has more than one feature missing, then the neighbors for that
    sample can be different depending on the particular feature being imputed.
    When the number of available neighbors is less than ``n_neighbors``, the
    training set average for that feature is used during imputation. When a
    sample has more than a ``feature_max_missing`` fraction of its features
    missing, then it is excluded from being a neighbor for imputation.

    .. versionadded:: 0.22

    Parameters
    ----------
    missing_values : number, string, np.nan, None (default=np.nan)
        The placeholder for the missing values. All occurrences of
        `missing_values` will be imputed.

    n_neighbors : int, optional (default=5)
        Number of neighboring samples to use for imputation.

    weights : str or callable, optional (default="uniform")
        Weight function used in prediction.  Possible values:

        - 'uniform' : uniform weights. All points in each neighborhood are
          weighted equally.
        - 'distance' : weight points by the inverse of their distance.
          in this case, closer neighbors of a query point will have a
          greater influence than neighbors which are further away.
        - callable : a user-defined function which accepts an
          array of distances, and returns an array of the same shape
          containing the weights.

    metric : str or callable, optional (default="nan_euclidean")
        Distance metric for searching neighbors. Possible values:
        - 'nan_euclidean'
        - callable : a user-defined function which conforms to the
          definition of _pairwise_callable(X, Y, metric, **kwds). The function
          accepts two arrays, X and Y, and a ``missing_values`` keyword in
          ``kwds`` and returns a scalar distance value.

    copy : boolean, optional(default=True)
        If True, a copy of X will be created. If False, imputation will
        be done in-place whenever possible.

    Attributes
    ----------
    statistics_ : array, shape (n_features, )
        The 1-D array contains the mean of each feature calculated using
        observed (i.e. non-missing) values. This is used for imputing
        missing values in samples that are either excluded from nearest
        neighbors search because they have too many missing features or
        all of the sample's k-nearest neighbors (i.e., the potential donors)
        also have the relevant feature value missing.

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
    >>> imputer = KNNImputer(n_neighbors=2)
    >>> imputer.fit_transform(X)
    array([[1. , 2. , 4. ],
           [3. , 4. , 3. ],
           [5.5, 6. , 5. ],
           [8. , 8. , 7. ]])
    """

    def __init__(self, missing_values=np.nan, n_neighbors=5,
                 weights="uniform", metric="nan_euclidean", copy=True):

        self.missing_values = missing_values
        self.n_neighbors = n_neighbors
        self.weights = weights
        self.metric = metric
        self.copy = copy

    def _calc_impute(self, dist, receivers_idx, X_col, fit_X_col,
                     non_missing_fix_X_col, statistic):
        """Helper function to calculate the imputation value of a single
        column.

        Parameters
        ----------
        dist : array-like, shape = (n_receivers, n_train_samples, )
            distance matrix between the receivers and the train samples

        receivers_idx : array-like, shape = (n_receivers, )
            indices of X_col to be imputed

        X_col : array-like, shape = (n_samples, )
            column to be imputed

        fit_X_col : array-like, shape = (n_train_samples, )
            column from training

        non_missing_fix_X_col : array-like, shape = (n_train_samples, )
            mask for missing values corresponding to ``fit_X_col``

        statistic : float
            statistic for column generated from fitting ``fit_X_col``
        """
        # Row index for receivers and potential donors (pdonors)
        potential_donors_idx = np.flatnonzero(non_missing_fix_X_col)

        # Use statistics if not enough donors are available
        if len(potential_donors_idx) < self.n_neighbors:
            return statistic

        # Get distance from potential donors
        dist_potential_donors = dist[:, potential_donors_idx]

        # Get donors
        donors_idx = np.argpartition(dist_potential_donors,
                                     self.n_neighbors - 1,
                                     axis=1)[:, :self.n_neighbors]

        # Get weight matrix from from distance matrix
        donors_dist = dist_potential_donors[
            np.arange(donors_idx.shape[0])[:, None], donors_idx]
        weight_matrix = _get_weights(donors_dist, self.weights)

        # Retrieve donor values and calculate kNN score
        donors = fit_X_col[potential_donors_idx].take(donors_idx)
        donors = np.ma.array(donors,
                             mask=_get_missing_mask(donors,
                                                    self.missing_values))

        return np.ma.average(donors, axis=1, weights=weight_matrix).data

    def fit(self, X, y=None):
        """Fit the imputer on X.

        Parameters
        ----------
        X : array-like, shape = (n_samples, n_features)
            Input data, where ``n_samples`` is the number of samples and
            ``n_features`` is the number of features.

        Returns
        -------
        self : object
        """

        # Check data integrity and calling arguments
        if not is_scalar_nan(self.missing_values):
            force_all_finite = True
        else:
            force_all_finite = "allow-nan"
            if self.metric not in _NAN_METRICS and not callable(self.metric):
                raise ValueError(
                    "The selected metric does not support NaN values")
        if self.n_neighbors <= 0:
            raise ValueError(
                "Expected n_neighbors > 0. Got {}".format(self.n_neighbors))

        X = check_array(X, accept_sparse=False, dtype=FLOAT_DTYPES,
                        force_all_finite=force_all_finite, copy=self.copy)

        mask = _get_missing_mask(X, self.missing_values)
        self.weights = _check_weights(self.weights)
        self._fit_X = X

        masked_X = np.ma.masked_array(X, mask=mask)
        mean_masked = masked_X.mean(axis=0)
        self.statistics_ = mean_masked.data
        self.statistics_[mean_masked.mask] = np.nan

        return self

    def transform(self, X):
        """Impute all missing values in X.

        Parameters
        ----------
        X : array-like, shape = (n_samples, n_features)
            The input data to complete.

        Returns
        -------
        X : array-like, shape = (n_samples, n_features)
            The imputed dataset.
        """

        check_is_fitted(self, ["statistics_"])
        if not is_scalar_nan(self.missing_values):
            force_all_finite = True
        else:
            force_all_finite = "allow-nan"
        X = check_array(X, accept_sparse=False, dtype=FLOAT_DTYPES,
                        force_all_finite=force_all_finite, copy=self.copy)

        if X.shape[1] != self._fit_X.shape[1]:
            raise ValueError("Incompatible dimension between the fitted "
                             "dataset and the one to be transformed")

        mask = _get_missing_mask(X, self.missing_values)
        if not np.any(mask):
            return X

        row_missing_idx = np.flatnonzero(mask.any(axis=1))

        # Pairwise distances between receivers and fitted samples
        dist = pairwise_distances(X[row_missing_idx, :], self._fit_X,
                                  metric=self.metric,
                                  squared=False,
                                  missing_values=self.missing_values,
                                  force_all_finite=False)

        # Maps from indices from X to indices in dist matrix
        dist_idx_map = np.zeros(X.shape[0], dtype=np.int)
        dist_idx_map[row_missing_idx] = np.arange(0, row_missing_idx.shape[0])

        invalid_mask = _get_missing_mask(self.statistics_, np.nan)
        valid_mask = np.logical_not(invalid_mask)
        valid_statistics_indexes = np.flatnonzero(valid_mask)

        mask_fit_X = _get_missing_mask(self._fit_X, self.missing_values)
        non_missing_fix_X = np.logical_not(mask_fit_X)
        # Find and impute missing
        for col in range(X.shape[1]):

            # column has no missing values
            if not np.any(mask[:, col]):
                continue

            # column was all missing during training
            if col not in valid_statistics_indexes:
                continue

            receivers_idx = np.where(mask[:, col])[0]
            dist_subset = dist[dist_idx_map[receivers_idx]]
            value = self._calc_impute(dist_subset, receivers_idx, X[:, col],
                                      self._fit_X[:, col],
                                      non_missing_fix_X[:, col],
                                      self.statistics_[col])
            X[receivers_idx, col] = value

        X = X[:, valid_statistics_indexes]
        return X

    def _more_tags(self):
        return {'allow_nan': True}
