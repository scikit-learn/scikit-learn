"""
Utility for making estimators frozen / un-trainable.
"""
# Author: Joel Nothman
# License: BSD

from .base import BaseEstimator, MetaEstimatorMixin
from .utils.metaestimators import if_delegate_has_method

__all__ = ['FreezeWrap']


class FreezeWrap(BaseEstimator, MetaEstimatorMixin):
    """Disable fitting and cloning for the wrapped estimator

    Wrapping an estimator in this freezes it, such that:

    * ``clone(FreezeWrap(estimator))`` will return the same model without
      clearing it
    * ``FreezeWrap(estimator).fit(...)`` will not call ``estimator.fit()``
    * ``FreezeWrap(estimator).fit_transform(X, y)`` will just return
      ``estimator.transform(X)``

    Read more in the :ref:`User Guide <freeze>`.

    Parameters
    ----------
    estimator : estimator

    Notes
    -----
    Any keyword arguments passed to ``fit_transform``, will *not*
    be passed on to ``transform`` (and similar for ``fit_predict``).
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
