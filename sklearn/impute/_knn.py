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


class KNNImputer(BaseEstimator, TransformerMixin):
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

    Parameters
    ----------
    missing_values : number, string, np.nan, None (default=np.nan)
        The placeholder for the missing values. All occurrences of
        `missing_values` will be imputed.

    n_neighbors : int, optional (default = 5)
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

    feature_max_missing : float, optional (default=0.5)
        The maximum fraction of features that can be missing
        before the sample is excluded from nearest neighbor imputation. These
        samples will not be considered a potential donor in ``fit``, and in
        ``transform`` their missing feature values will be imputed to be the
        feature mean for the entire dataset.

    sample_max_missing : float, optional (default=0.8)
        The maximum fraction of samples that can be missing for any feature
        beyond which an error is raised.

    copy : boolean, optional(default=True)
        If True, a copy of X will be created. If False, imputation will
        be done in-place whenever possible. Note that, if metric is
        "nan_euclidean" and copy=False then missing_values in the
        input matrix X will be overwritten with zeros.

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
    >>> imputer = KNNImputer(n_neighbors=2, weights="uniform")
    >>> imputer.fit_transform(X)
    array([[1. , 2. , 4. ],
           [3. , 4. , 3. ],
           [5.5, 6. , 5. ],
           [8. , 8. , 7. ]])
    """

    def __init__(self, missing_values=np.nan, n_neighbors=5,
                 weights="uniform", metric="nan_euclidean",
                 feature_max_missing=0.5, sample_max_missing=0.8, copy=True):

        self.missing_values = missing_values
        self.n_neighbors = n_neighbors
        self.weights = weights
        self.metric = metric
        self.feature_max_missing = feature_max_missing
        self.sample_max_missing = sample_max_missing
        self.copy = copy

    def _impute(self, dist, receivers_idx,
                X_col, fit_X_col, mask_fit_X_col, statistic):
        """Helper function to impute the missing values of a single column.

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

        mask_fit_X_col : array-like, shape = (n_train_samples, )
            mask for missing values corresponding to ``fit_X_col``

        statistic : float
            statistic for column generated from fitting ``fit_X_col``
        """
        if receivers_idx.shape[0] == 0:
            return

        # Row index for receivers and potential donors (pdonors)
        potential_donors_idx = np.where(~mask_fit_X_col)[0]

        # Impute using column mean if n_neighbors are not available
        if len(potential_donors_idx) < self.n_neighbors:
            X_col[receivers_idx] = statistic
            return

        # Get distance from potential donors
        dist_potential_donors = dist[:, potential_donors_idx]

        # Get donors
        donors_idx = np.argpartition(dist_potential_donors,
                                     self.n_neighbors - 1,
                                     axis=1)[:, :self.n_neighbors]

        # Get weights or None
        donors_dist = dist_potential_donors[
            np.arange(donors_idx.shape[0])[:, None], donors_idx]
        weight_matrix = _get_weights(donors_dist, self.weights)

        # Retrieve donor values and calculate kNN score
        donors = fit_X_col[potential_donors_idx].take(donors_idx)
        donors = np.ma.array(donors,
                             mask=_get_missing_mask(donors,
                                                    self.missing_values))

        # Final imputation
        imputed = np.ma.average(donors, axis=1, weights=weight_matrix)
        X_col[receivers_idx] = imputed.data

    def fit(self, X, y=None):
        """Fit the imputer on X.

        Parameters
        ----------
        X : array-like, shape (n_samples, n_features, )
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
        X = check_array(X, accept_sparse=False, dtype=FLOAT_DTYPES,
                        force_all_finite=force_all_finite, copy=self.copy)

        # Check if sufficient neighboring samples available
        if X.shape[0] < self.n_neighbors:
            raise ValueError("There are only {} samples, but n_neighbors={}"
                             .format(X.shape[0], self.n_neighbors))

        # Check if % missing in any samples > sample_max_missing
        mask = _get_missing_mask(X, self.missing_values)
        if np.any(mask.sum(axis=0) > (X.shape[0] * self.sample_max_missing)):
            raise ValueError("Some columns have more than {}% missing values"
                             .format(self.sample_max_missing * 100))

        self.weights = _check_weights(self.weights)
        self.fit_X_ = X[
            mask.sum(axis=1) <= (mask.shape[1] * self.feature_max_missing)]
        self.statistics_ = np.ma.array(X, mask=mask).mean(axis=0).data

        return self

    def transform(self, X):
        """Impute all missing values in X.

        Parameters
        ----------
        X : array-like, shape = [n_samples, n_features]
            The input data to complete.

        Returns
        -------
        X : array-like, shape = [n_samples, n_features]
            The imputed dataset.
        """

        check_is_fitted(self, ["fit_X_", "statistics_"])
        if not is_scalar_nan(self.missing_values):
            force_all_finite = True
        else:
            force_all_finite = "allow-nan"
        X = check_array(X, accept_sparse=False, dtype=FLOAT_DTYPES,
                        force_all_finite=force_all_finite, copy=self.copy)

        if X.shape[1] != self.fit_X_.shape[1]:
            raise ValueError("Incompatible dimension between the fitted "
                             "dataset and the one to be transformed")

        mask = _get_missing_mask(X, self.missing_values)
        orig_X_shape = X.shape
        if not np.any(mask):
            return X

        # Check for excessive missingness in rows
        row_total_missing = mask.sum(axis=1)
        bad_rows = row_total_missing > (mask.shape[1] *
                                        self.feature_max_missing)
        if np.any(bad_rows):
            X_bad = X[bad_rows]
            X = X[~bad_rows]
            mask = mask[~bad_rows]
            row_total_missing = mask.sum(axis=1)

        row_misses = row_total_missing.astype(np.bool)
        if np.any(row_misses):
            # Pairwise distances between receivers and fitted samples
            dist = pairwise_distances(X[row_misses], self.fit_X_,
                                      metric=self.metric,
                                      squared=False,
                                      missing_values=self.missing_values)

            # Maps from indices from X to indices in dist matrix
            dist_idx_map = np.zeros(X.shape[0], dtype=np.int)
            dist_idx_map[row_misses] = np.arange(0, row_misses.sum(0))

            # Find and impute missing
            mask_fit_X = _get_missing_mask(self.fit_X_, self.missing_values)
            for col in range(X.shape[1]):

                if not np.any(mask[:, col]):
                    continue

                receivers_idx = np.where(mask[:, col])[0]
                dist_subset = dist[dist_idx_map[receivers_idx]]
                self._impute(dist_subset, receivers_idx, X[:, col],
                             self.fit_X_[:, col], mask_fit_X[:, col],
                             self.statistics_[col])

        # Merge deficient rows and mean impute their missing values
        if np.any(bad_rows):
            missing_mask = _get_missing_mask(X_bad, self.missing_values)
            bad_missing_index = np.where(missing_mask)
            X_bad[bad_missing_index] = self.statistics_[bad_missing_index[1]]
            X_merged = np.empty(orig_X_shape)
            X_merged[bad_rows] = X_bad
            X_merged[~bad_rows] = X
            X = X_merged
        return X

    def _more_tags(self):
        return {'allow_nan': True}
