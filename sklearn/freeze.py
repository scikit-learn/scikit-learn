"""
Utility for making estimators frozen / un-trainable.
"""

from .base import BaseEstimator
from .utils.metaestimators import if_delegate_has_method


class FreezeWrap(BaseEstimator):
    """Disables fitting and cloning (clearing model) for the wrapped estimator

    Parameters
    ----------
    estimator : estimator
    """

    def __init__(self, estimator):
        self.estimator = estimator

    def fit(self, X, y=None, **kwargs):
        """Return self

        Parameters
        ----------
        X
            ignored
        y : optional
            ignored
        kwargs : optional
            ignored
        """
        return self

    @if_delegate_has_method(delegate='estimator')
    def fit_transform(self, X, y=None, **kwargs):
        """Execute transform on estimator

        Parameters
        ----------
        X
            data to transform
        y : optional
            ignored
        kwargs : ignored
            ignored
        """
        return self.estimator.transform(X)

    @if_delegate_has_method(delegate='estimator')
    def fit_predict(self, X, y=None, **kwargs):
        """Execute predict on estimator

        Parameters
        ----------
        X
            data to predict
        y : optional
            ignored
        kwargs : ignored
            ignored
        """
        return self.estimator.predict(X)

    @if_delegate_has_method(delegate='estimator')
    def transform(self, *args, **kwargs):
        """Execute estimator's equivalent method
        """
        return self.estimator.transform(*args, **kwargs)

    @if_delegate_has_method(delegate='estimator')
    def decision_function(self, *args, **kwargs):
        """Execute estimator's equivalent method
        """
        return self.estimator.decision_function(*args, **kwargs)

    @if_delegate_has_method(delegate='estimator')
    def predict(self, *args, **kwargs):
        """Execute estimator's equivalent method
        """
        return self.estimator.predict(*args, **kwargs)

    @if_delegate_has_method(delegate='estimator')
    def predict_log_proba(self, *args, **kwargs):
        """Execute estimator's equivalent method
        """
        return self.estimator.predict_log_proba(*args, **kwargs)

    @if_delegate_has_method(delegate='estimator')
    def predict_proba(self, *args, **kwargs):
        """Execute estimator's equivalent method
        """
        return self.estimator.predict_proba(*args, **kwargs)

    @property
    def _estimator_type(self):
        return self.estimator._estimator_type

    @property
    def classes_(self):
        return self.estimator.classes_

    @property
    def _pairwise(self):
        # check if first estimator expects pairwise input
        return getattr(self.estimator, '_pairwise', False)
