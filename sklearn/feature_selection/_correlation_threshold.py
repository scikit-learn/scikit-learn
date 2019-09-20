import numpy as np
from scipy.sparse import issparse

from ..base import BaseEstimator
from .base import SelectorMixin
from ..utils import check_array
from ..utils.validation import check_is_fitted
from ..utils.extmath import safe_sparse_dot
from ..utils.sparsefuncs import min_max_axis, mean_variance_axis


class CorrelationThreshold(BaseEstimator, SelectorMixin):
    """Feature selector that removes features with high correlation.

    This unsupervised feature selection algorithm considers feature pairs with
    correlation higher than the `threshold` and removes the feature that has
    the highest mean correlation with the other features. This process is
    continued till there are no correlations higher than `threshold`.

    Read more in the :ref:`User Guide <correlation_threshold>`.

    Parameters
    ----------
    threshold : float, default=0.9
        Features with a training-set correlation higher than this threshold
        will be removed.

    algorithm : {'pearson', 'spearmanr'}, default='pearson'
        Correlation type to use. Dense arrays can use any type. 'pearson' is
        the only type that supports sparse matrices.

    Attributes
    ----------
    support_mask_ : bool ndarray of shape (n_features,)
        Boolean mask for features to keep.

    Examples
    --------
    The following example shows how one of the first highly correlated features
    are removed, while the uncorrelated features remains:

    >>> import numpy as np
    >>> from sklearn.feature_selection import CorrelationThreshold
    >>> X = np.array([[0.0, 1.0, 2.0], [1.1, 2.0, 3.0], [0.5, 10.1, 1.1]]).T
    >>> selector = CorrelationThreshold(threshold=0.9)
    >>> selector.fit_transform(X)
    array([[ 1.1,  0.5],
           [ 2. , 10.1],
           [ 3. ,  1.1]])

    References
    ----------
    Max Kuhn and Kjell Johnson, "Applied Predictive Modeling", Springer, 2013
    (Section 3.5)
    """
    def __init__(self, threshold=0.9, algorithm='pearson'):
        self.threshold = threshold
        self.algorithm = algorithm

    def fit(self, X, y=None):
        """Learn empirical variances from X.

        Parameters
        ----------
        X : {array-like, sparse matrix} of shape (n_samples, n_features)
            Training set to compute correlations.

        y : ignored
            Not used, present here for API consistency by convention.

        Returns
        -------
        self : `CorrelationThreshold`
        """
        if not (0.0 <= self.threshold <= 1.0):
            raise ValueError("threshold must be in [0.0, 1.0], got {}".format(
                             self.threshold))
        if issparse(X) and self.algorithm != 'pearson':
            raise ValueError("only pearson correlation is supported with "
                             "sparse matrices")

        X = check_array(X, accept_sparse=['csc', 'csr'], dtype=[np.float64,
                                                                np.float32])

        n_features = X.shape[1]
        if self.threshold == 1.0 or (1 in X.shape):
            self.support_mask_ = np.ones(n_features, dtype=np.bool)
            return self

        # get constant features
        if issparse(X):
            mins, maxes = min_max_axis(X, axis=0)
            peak_to_peaks = maxes - mins
            constant_mask = np.isclose(peak_to_peaks, 0.0)

            # sparse correlation
            mu, sparse_var = mean_variance_axis(X, 0)
            X_corr = _sparse_correlation(X, mu, ~constant_mask)
        else:
            peak_to_peaks = np.ptp(X, axis=0)
            constant_mask = np.isclose(peak_to_peaks, 0.0)
            X_corr = np.corrcoef(X, rowvar=False)

        np.fabs(X_corr, out=X_corr)

        # Removes constant features from support_mask
        support_mask = np.ones(n_features, dtype=np.bool)
        upper_idx = np.triu_indices(n_features, 1)

        non_constant_features = n_features
        for i in np.flatnonzero(constant_mask):
            feat_remove_mask = np.logical_and(upper_idx[0] != i,
                                              upper_idx[1] != i)
            upper_idx = (upper_idx[0][feat_remove_mask],
                         upper_idx[1][feat_remove_mask])
            support_mask[i] = False
            non_constant_features -= 1

        for _ in range(non_constant_features - 1):
            max_idx = np.argmax(X_corr[upper_idx])
            feat1, feat2 = upper_idx[0][max_idx], upper_idx[1][max_idx]
            cur_corr = X_corr[feat1, feat2]

            # max correlation is lower than threshold
            if cur_corr < self.threshold:
                break

            # Temporary remove both features to calculate the mean with other
            # features. One of the features will be selected.
            support_mask[[feat1, feat2]] = False

            # if there are no other features to compare, keep the feature with
            # the most variance
            if np.all(~support_mask):
                if issparse(X):
                    var = sparse_var[[feat1, feat2]]
                else:
                    var = np.var(X[:, [feat1, feat2]], axis=0)

                if var[0] > var[1]:
                    support_mask[feat1] = True
                else:
                    support_mask[feat2] = True
                break

            # means with other features
            feat1_mean = np.mean(X_corr[feat1, support_mask])
            feat2_mean = np.mean(X_corr[feat2, support_mask])

            # feature with lower mean is kept
            if feat1_mean > feat2_mean:
                support_mask[feat2] = True
                feat_to_remove = feat1
            else:
                support_mask[feat1] = True
                feat_to_remove = feat2

            # Removes the removed feature from consideration
            upper_idx_to_keep = np.logical_and(upper_idx[0] != feat_to_remove,
                                               upper_idx[1] != feat_to_remove)
            upper_idx = (upper_idx[0][upper_idx_to_keep],
                         upper_idx[1][upper_idx_to_keep])

        self.support_mask_ = support_mask
        return self

    def _get_support_mask(self):
        check_is_fitted(self)
        return self.support_mask_


def _sparse_correlation(X, mu, non_constant_mask):
    """Calcuate Pearson correlation for sparse matrices

    Parameters
    ----------
    X : sparse matrix of shape (n_samples, n_features)
        Matrix to find correlation on.

    mu : ndarray of shape (n_features,)
        Mean of feature columns.

    non_constant_mask : ndarray of shape (n_features,)
        Boolean mask for non constant features.

    Returns
    -------
    correlation matrix : ndarray of shape (n_features, n_features)
    """
    X_diff = X - mu[None, :]
    X_corr = safe_sparse_dot(X_diff.T, X_diff, dense_output=True)
    stddev = np.sqrt(np.diag(X_corr))

    X_corr[non_constant_mask, :] /= stddev[non_constant_mask][:, None]
    X_corr[:, non_constant_mask] /= stddev[non_constant_mask][None, :]
    return X_corr
