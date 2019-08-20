import numpy as np
from scipy.sparse import issparse

from ..base import BaseEstimator
from .base import SelectorMixin
from ..utils import check_array
from ..utils.validation import check_is_fitted
from ..utils.extmath import safe_sparse_dot
from ..utils.sparsefuncs import min_max_axis


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

    Attributes
    ----------
    support_mask_ : bool ndarray of shape (n_features, )
        Boolean mask for features to keep

    Examples
    --------
    The following example shows how one of the first highly correlated features
    are removed, while the uncorrelated feature remains:

    >>> import numpy as np
    >>> from sklearn.feature_selection import CorrelationThreshold
    >>> X = np.array([[0.0, 1.0, 2.0], [1.1, 2.0, 3.0], [0.5, 10.1, 1.1]]).T
    >>> selector = CorrelationThreshold()
    >>> selector.fit_transform(X)
    array([[ 1.1,  0.5],
           [ 2. , 10.1],
           [ 3. ,  1.1]])
    """
    def __init__(self, threshold=0.9):
        self.threshold = threshold

    def fit(self, X, y=None):
        """Learn empirical variances from X.

        Parameters
        ----------
        X : {array-like, sparse matrix} shape of (n_samples, n_features)
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

        X = check_array(X, accept_sparse=['csc', 'csr'], dtype=float)

        n_features = X.shape[1]
        if self.threshold == 1.0:
            self.support_mask_ = np.ones(n_features, dtype=np.bool)
            return self

        support_mask = np.ones(n_features, dtype=np.bool)
        upper_idx = np.triu_indices(n_features, 1)

        # get constant features
        if issparse(X):
            mins, maxes = min_max_axis(X, axis=0)
            peak_to_peaks = maxes - mins
        else:
            peak_to_peaks = np.ptp(X, axis=0)
        ptp_mask = peak_to_peaks != 0

        X_corr = safe_sparse_dot(X.T, X, dense_output=True)
        # only divide when feature is non constant
        stddev = np.sqrt(np.diag(X_corr))[ptp_mask]
        X_corr[ptp_mask, :] /= stddev[:, None]
        X_corr[:, ptp_mask] /= stddev[None, :]

        np.fabs(X_corr, out=X_corr)

        # Removes constant features from support_mask
        non_constant_features = n_features
        for i, ptp in enumerate(peak_to_peaks):
            if ptp == 0.0:
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

            support_mask[[feat1, feat2]] = False

            # if there are no other features to compare to then keep the keep
            # the first feature
            if np.all(~support_mask):
                support_mask[feat1] = True
                break

            # means with other features
            feat1_mean = np.mean(X_corr[feat1, support_mask])
            feat2_mean = np.mean(X_corr[feat2, support_mask])

            # feature with higher mean is removed
            if feat1_mean > feat2_mean:
                feat_to_keep = feat2
                feat_to_remove = feat1
            else:
                feat_to_keep = feat1
                feat_to_remove = feat2

            support_mask[feat_to_keep] = True

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
