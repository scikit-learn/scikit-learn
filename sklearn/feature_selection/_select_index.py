# Author: Lars Buitinck
# License: 3-clause BSD

import numpy as np
from ..base import BaseEstimator
from ._base import SelectorMixin
from ..utils.validation import check_is_fitted


class SelectColumnsByIndicies(SelectorMixin, BaseEstimator):
    """Feature selector given columns indicies.

    This feature selection algorithm looks only at the features (X), not the
    desired outputs (y), and can thus be used for unsupervised learning.

    Parameters
    ----------
    indicies : array
        Indicies of columns/features to be preserved.

    Attributes
    ----------
    n_features_in_ : int
        Number of features seen during :term:`fit`.

        .. versionadded:: 0.24

    feature_names_in_ : ndarray of shape (`n_features_in_`,)
        Names of features seen during :term:`fit`. Defined only when `X`
        has feature names that are all strings.

        .. versionadded:: 1.0

    See Also
    --------
    SelectFromModel: Meta-transformer for selecting features based on
        importance weights.
    SelectPercentile : Select features according to a percentile of the highest
        scores.
    SequentialFeatureSelector : Transformer that performs Sequential Feature
        Selection.

    Notes
    -----
    Allows NaN in the input.

    Examples
    --------
    The following dataset has integer features, two of which are the same
    in every sample. These are removed with the default setting for threshold::

        >>> X = [[0, 2, 0, 3], [0, 1, 4, 3], [0, 1, 1, 3]]
        >>> selector = SelectColumnsByIndicies(indicies=[1, 2])
        >>> selector.fit_transform(X)
        array([[2, 0],
               [1, 4],
               [1, 1]])
    """

    def __init__(self, indicies):
        self.indicies = np.array(indicies)

    def fit(self, X, y=None):
        """Select columns given indicies.

        Parameters
        ----------
        X : {array-like, sparse matrix}, shape (n_samples, n_features)
            Data from which select columns given indicies.

        y : any, default=None
            Ignored. This parameter exists only for compatibility with
            sklearn.pipeline.Pipeline.

        Returns
        -------
        self : object
            Returns the instance itself.
        """

        X = self._validate_data(
            X,
            accept_sparse=("csr", "csc"),
            dtype=np.float64,
            force_all_finite="allow-nan",
        )

        support_mask = np.zeros(self.n_features_in_, dtype=bool)
        support_mask[self.indicies] = True

        self.support_mask = support_mask

        return self

    def _get_support_mask(self):
        check_is_fitted(self)

        return self.support_mask

    def _more_tags(self):
        return {"allow_nan": True}
