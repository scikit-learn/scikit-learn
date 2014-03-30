# Authors: Gilles Louppe, Mathieu Blondel, Maheshakya Wijewardena
# License: BSD 3 clause

import numpy as np

from ..base import TransformerMixin, BaseEstimator, clone
from ..externals import six
from ..utils import safe_mask, atleast2d_or_csc, deprecated
from .base import SelectorMixin


class _LearntSelectorMixin(TransformerMixin):
    # Note because of the extra threshold parameter in transform, this does
    # not naturally extend from SelectorMixin

    """Transformer mixin selecting features based on importance weights.

    This implementation can be mixin on any estimator that exposes a
    ``feature_importances_`` or ``coef_`` attribute to evaluate the relative
    importance of individual features for feature selection.
    """
    @deprecated('Support to use directly on estimators '
                'will be removed in version 0.17')
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
        X = atleast2d_or_csc(X)

        if isinstance(self, SelectorMixin
                      ) and isinstance(self, TransformerMixin):
            importances = self._set_importances(self.estimator, X)
            threshold = self._set_threshold(self.estimator, self.threshold)
        else:
            importances = self._set_importances(self.estimator, X)
            threshold = self._set_threshold(self, threshold)

        if isinstance(threshold, six.string_types):
            if "*" in threshold:
                scale, reference = threshold.split("*")
                scale = float(scale.strip())
                reference = reference.strip()

                if reference == "median":
                    reference = np.median(importances)
                elif reference == "mean":
                    reference = np.mean(importances)
                else:
                    raise ValueError("Unknown reference: " + reference)

                threshold = scale * reference

            elif threshold == "median":
                threshold = np.median(importances)

            elif threshold == "mean":
                threshold = np.mean(importances)

        else:
            threshold = float(threshold)

        self.threshold_ = threshold

        # Selection
        try:
            mask = importances >= threshold
        except TypeError:
            # Fails in Python 3.x when threshold is str;
            # result is array of True
            raise ValueError("Invalid threshold: all features are discarded.")

        if np.any(mask):
            mask = safe_mask(X, mask)
            self._support_mask = mask
            return X[:, mask]
        else:
            raise ValueError("Invalid threshold: all features are discarded.")

    @staticmethod
    def _set_importances(estimator, X):
        # Retrieve importance vector
        if hasattr(estimator, "feature_importances_"):
            importances = estimator.feature_importances_
            if importances is None:
                raise ValueError("Importance weights not computed. Please set"
                                 " the compute_importances parameter before "
                                 "fit.")

        elif hasattr(estimator, "coef_"):
            if estimator.coef_.ndim == 1:
                importances = np.abs(estimator.coef_)

            else:
                importances = np.sum(np.abs(estimator.coef_), axis=0)

        else:
            raise ValueError("Missing `feature_importances_` or `coef_`"
                             " attribute, did you forget to set the "
                             "estimator's parameter to compute it?")
        if len(importances) != X.shape[1]:
            raise ValueError("X has different number of features than"
                             " during model fitting.")

        return importances

    @staticmethod
    def _set_threshold(estimator, threshold):
        # Retrieve threshold
        if threshold is None:
            if hasattr(estimator, "penalty") and estimator.penalty == "l1":
                # the natural default threshold is 0 when l1 penalty was used
                threshold = getattr(estimator, "threshold", 1e-5)
            else:
                threshold = getattr(estimator, "threshold", "mean")

        return threshold


class SelectFromModel(BaseEstimator, _LearntSelectorMixin,
                      SelectorMixin):

    """Meta-transformer for selecting features based on importance
       weights.

    Parameters
    ----------
    estimator : object or None(default=None)
        The base estimator from which the transformer is built.
        If None, then a value error is raised.

    threshold : string, float or None, optional (default=None)
            The threshold value to use for feature selection. Features whose
            importance is greater or equal are kept while the others are
            discarded. If "median" (resp. "mean"), then the threshold value is
            the median (resp. the mean) of the feature importances. A scaling
            factor (e.g., "1.25*mean") may also be used. If None and if
            available, the object attribute ``threshold`` is used. Otherwise,
            "mean" is used by default.

    Attributes
    ----------
    `estimator_`: an estimator
        The base estimator from which the transformer is built.

    `threshold_`: float
        The threshold value used for feature selection.
    """

    def __init__(self, estimator, threshold=None):
        self.estimator = estimator
        self.threshold = threshold

    def _get_support_mask(self):
        return self._support_mask

    def fit(self, X, y, **fit_params):
        """Trains the estimator in order to perform transformation
           set (X, y).

        Parameters
        ----------
        X : array-like of shape = [n_samples, n_features]
            The training input samples.

        y : array-like, shape = [n_samples]
            The target values (integers that correspond to classes in
            classification, real numbers in regression).

        **fit_params : Other estimator specific parameters

        Returns
        -------
        self : object
            Returns self.
        """

        # Validate and make estimator
        if self.estimator is not None:
            self.estimator_ = self.estimator
        else:
            raise ValueError("estimator cannot be None")

        self.estimator = clone(self.estimator_)

        # Fit the estimator
        self.estimator.fit(X, y, **fit_params)

        return self
