# Author: Lars Buitinck
# License: 3-clause BSD

import numpy as np
from ..base import BaseEstimator
from .base import SelectorMixin
from ..utils import check_array
from ..utils.sparsefuncs import mean_variance_axis, min_max_axis
from ..utils.validation import check_is_fitted


class VarianceThreshold(SelectorMixin, BaseEstimator):
    """Feature selector that removes all low-variance features.

    This feature selection algorithm looks only at the features (X), not the
    desired outputs (y), and can thus be used for unsupervised learning.

    Read more in the :ref:`User Guide <variance_threshold>`.

    Parameters
    ----------
    threshold : float, optional
        Features with a training-set variance lower than this threshold will
        be removed. The default is to keep all features with non-zero variance,
        i.e. remove the features that have the same value in all samples.

    Attributes
    ----------
    variances_ : array, shape (n_features,)
        Variances of individual features.

    Examples
    --------
    The following dataset has integer features, two of which are the same
    in every sample. These are removed with the default setting for threshold::

        >>> X = [[0, 2, 0, 3], [0, 1, 4, 3], [0, 1, 1, 3]]
        >>> selector = VarianceThreshold()
        >>> selector.fit_transform(X)
        array([[2, 0],
               [1, 4],
               [1, 1]])
    """

    def __init__(self, threshold=0.):
        self.threshold = threshold

    def fit(self, X, y=None):
        """Learn empirical variances from X.

        Parameters
        ----------
        X : {array-like, sparse matrix}, shape (n_samples, n_features)
            Sample vectors from which to compute variances.

        y : any
            Ignored. This parameter exists only for compatibility with
            sklearn.pipeline.Pipeline.

        Returns
        -------
        self
        """
        X = check_array(X, ('csr', 'csc'), dtype=np.float64)

        if hasattr(X, "toarray"):   # sparse matrix
            _, self.variances_ = mean_variance_axis(X, axis=0)
            if self.threshold == 0:
                mins, maxes = min_max_axis(X, axis=0)
                peak_to_peaks = maxes - mins
        else:
            self.variances_ = np.var(X, axis=0)
            if self.threshold == 0:
                peak_to_peaks = np.ptp(X, axis=0)

        if self.threshold == 0:
            # Use peak-to-peak to avoid numeric precision issues
            # for constant features
            self.variances_ = np.minimum(self.variances_, peak_to_peaks)

        if np.all(self.variances_ <= self.threshold):
            msg = "No feature in X meets the variance threshold {0:.5f}"
            if X.shape[0] == 1:
                msg += " (X contains only one sample)"
            raise ValueError(msg.format(self.threshold))

        return self

    def _get_support_mask(self):
        check_is_fitted(self)

        return self.variances_ > self.threshold
