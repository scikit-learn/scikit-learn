# Authors: Gilles Louppe, Mathieu Blondel
# License: BSD 3 clause

import numpy as np

from ..base import TransformerMixin, BaseEstimator

from .etc import SelectByScoreMixin


class _Selector(BaseEstimator, SelectByScoreMixin):

    def __init__(self, estimator, minimum=None, maximum=None, scaling=None,
                 limit=None):
        self.estimator = estimator
        self.minimum = minimum
        self.maximum = maximum
        self.limit = limit
        self.scaling = scaling
        try:
            self._fit()
        except ValueError:
            # assume underlying estimator needs fitting
            pass

    def fit(self, X, y=None):
        self.estimator.fit(X, y)
        self._fit(X)
        return self

    def _fit(self, X=None):
        super(_Selector, self)._fit(self._get_feature_importances(), X)

    def _get_support_mask(self, minimum=None, maximum=None, scaling=None,
                          limit=None):

        return super(_Selector, self)._get_support_mask(
            minimum=minimum, maximum=maximum, scaling=scaling, limit=limit)


class SelectorMixin(TransformerMixin):
    """Transformer mixin selecting features based on importance weights.

    This implementation can be mixin on any estimator that exposes a
    ``feature_importances_`` or ``coef_`` attribute to evaluate the relative
    importance of individual features for feature selection.
    """

    def make_transformer(self):
        if getattr(self, 'penalty', None) == 'l1':
            # the natural default threshold is 0 when l1 penalty was used
            default_threshold = 1e-5
        else:
            default_threshold = 'mean'
        return _Selector(self, minimum=default_threshold)

    def _get_feature_importances(self, X=None, y=None):
        if hasattr(self, "feature_importances_"):
            importances = self.feature_importances_
            if importances is None:
                raise ValueError("Importance weights not computed. Please set"
                                 " the compute_importances parameter before "
                                 "fit.")

        elif hasattr(self, "coef_"):
            if self.coef_.ndim == 1:
                importances = np.abs(self.coef_)

            else:
                importances = np.sum(np.abs(self.coef_), axis=0)

        else:
            raise ValueError("Missing `feature_importances_` or `coef_`"
                             " attribute, did you forget to set the "
                             "estimator's parameter to compute it?")
        return importances

    def transform(self, X, threshold=None):
        """Reduce X to its most important features.

        Parameters
        ----------
        X : array or scipy sparse matrix of shape [n_samples, n_features]
            The input samples.

        threshold : string, float or None, optional (default=None)
            The threshold value to use for feature selection. Features whose
            importance is greater or equal are kept while the others are
            discarded. If "median" (resp. "mean"), then the threshold value is
            the median (resp. the mean) of the feature importances. A scaling
            factor (e.g., "1.25*mean") may also be used. If None and if
            available, the object attribute ``threshold`` is used. Otherwise,
            "mean" is used by default.

        Returns
        -------
        X_r : array of shape [n_samples, n_selected_features]
            The input samples with only the selected features.
        """
        transformer = self.make_transformer()
        if threshold is not None:
            transformer.minimum = threshold
        return transformer.fit_transform(X)
