# Authors: Gilles Louppe, Mathieu Blondel, Maheshakya Wijewardena
# License: BSD 3 clause

import numpy as np

from ..base import (TransformerMixin, BaseEstimator, clone, SelectorMixin,
    MetaEstimatorMixin)
from ..externals import six

from ..utils import safe_mask, check_array, deprecated
from ..utils.validation import NotFittedError, check_is_fitted


def _get_feature_importances(estimator, X):
    """Retrieve or aggregate feature importances from estimator"""
    if hasattr(estimator, "feature_importances_"):
        importances = estimator.feature_importances_

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


def _calculate_threshold(estimator, importances, threshold):
    """Interpret the threshold value"""

    if threshold is None:
        # determine default from estimator
        if hasattr(estimator, "penalty") and estimator.penalty == "l1":
            # the natural default threshold is 0 when l1 penalty was used
            threshold = 1e-5
        else:
            threshold = "mean"

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

    return threshold


class _LearntSelectorMixin(TransformerMixin):
    # Note because of the extra threshold parameter in transform, this does
    # not naturally extend from SelectorMixin
    """Transformer mixin selecting features based on importance weights.

    This implementation can be mixin on any estimator that exposes a
    ``feature_importances_`` or ``coef_`` attribute to evaluate the relative
    importance of individual features for feature selection.
    """
    @deprecated('Support to use estimators as feature selectors will be '
                'removed in version 0.17. Use SelectFromModel instead.')
    def transform(self, X, threshold=None):
        """Reduce X to its most important features.

        Uses ``coef_`` or ``feature_importances_`` to determine the most
        important features.  For models with a ``coef_`` for each class, the
        absolute sum over the classes is used.

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
        check_is_fitted(self, ('coef_', 'feature_importances_'), 
                        all_or_any=any)

        X = check_array(X, 'csc')
        importances = _get_feature_importances(self, X)

        if threshold is None:
            threshold = getattr(self, 'threshold', None)
        threshold = _calculate_threshold(self, importances, threshold)

        # Selection
        try:
            mask = importances >= threshold
        except TypeError:
            # Fails in Python 3.x when threshold is str;
            # result is array of True
            raise ValueError("Invalid threshold: all features are discarded.")

        if np.any(mask):
            mask = safe_mask(X, mask)
            return X[:, mask]
        else:
            raise ValueError("Invalid threshold: all features are discarded.")

    @staticmethod
    def _set_importances(estimator, X):
        # Retrieve importance vector
        if hasattr(self, "feature_importances_"):
            importances = self.feature_importances_

        elif hasattr(self, "coef_"):
            if self.coef_ is None:
                msg = "This model is not fitted yet. Please call fit() first" 
                raise NotFittedError(msg)

            if self.coef_.ndim == 1:
                importances = np.abs(self.coef_)

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


class SelectFromModel(BaseEstimator, SelectorMixin):
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

    warm_start : bool, optional
        When set to True, reuse the solution of the previous call to fit as
        initialization, otherwise, just erase the previous solution.

    Attributes
    ----------
    `estimator_`: an estimator
        The base estimator from which the transformer is built.

    `scores_`: array, shape=(n_features,)
        The importance of each feature according to the fit model.

    `threshold_`: float
        The threshold value used for feature selection.
    """

    def __init__(self, estimator, threshold=None, warm_start=False):
        self.estimator = estimator
        self.threshold = threshold
        self.warm_start = warm_start

    def _get_support_mask(self):
        self.threshold_ = _calculate_threshold(self.estimator, self.scores_,
                                               self.threshold)
        return self.scores_ >= self.threshold_

    def fit(self, X, y, **fit_params):
        """Fit the SelectFromModel meta-transformer.

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
        if not (self.warm_start and hasattr(self,"estimator_")):
            self.estimator_ = clone(self.estimator)

        self.estimator_.fit(X, y, **fit_params)
        self.scores_ = _get_feature_importances(self.estimator_, X)
        return self

    def partial_fit(self, X, y, **fit_params):
        """Fit the SelectFromModel meta-transformer only once.

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
        if not hasattr(self, "estimator_"):
            self.estimator_ = clone(self.estimator)
        self.estimator_.partial_fit(X, y, **fit_params)
        self.scores_ = _get_feature_importances(self.estimator_, X)
        return self
