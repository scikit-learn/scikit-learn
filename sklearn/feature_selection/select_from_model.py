"""Meta-transformer."""

# Author: Maheshakya Wijewardena <maheshakya.10@cse.mrt.ac.lk>
# License: BSD 3 clause

from ..base import BaseEstimator, clone
from ..utils import check_arrays

from .from_model import _LearntSelectorMixin


class SelectFromModel(BaseEstimator, _LearntSelectorMixin):

    """Meta-transformer for selecting features based on importance
       weights.

    Parameters
    ----------
    estimator : object or None(default=None)
        The base estimator from which the transformer is built.
        If None, then a value error is raised.

    estimator_params : list of strings
        The list of attributes to use as parameters when instantiating a
        new base estimator. If none are given, default parameters are used.

    Attributes
    ----------
    `estimator_`: an estimators
        The base estimator from which the transformer is built.
    """

    def _validate_estimator(self, default=None):
        """Check the estimator and set the `estimator_` attribute."""
        if self.estimator is not None:
            self.estimator_ = self.estimator
        else:
            self.estimator_ = default

        if self.estimator_ is None:
            raise ValueError("estimator cannot be None")

    def _make_estimator(self):
        """Make and configure a copy of the `estimator_` attribute."""

        estimator = clone(self.estimator_)
        estimator.set_params(**dict((p, getattr(self, p))
                                    for p in self.estimator_params))

        return estimator

    def __init__(self, estimator=None, estimator_params=tuple()):
        self.estimator = estimator
        self.estimator_params = estimator_params

        self._validate_estimator()
        self.estimator = self._make_estimator()

        if hasattr(self.estimator, "penalty"):
            self.penalty = self.estimator.penalty

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

        # Convert data
        X, y = check_arrays(X, y)

        # Fit the estimator
        self.estimator.fit(X, y, **fit_params)

        if hasattr(self.estimator, "feature_importances_"):
            self.feature_importances_ = self.estimator.feature_importances_

        elif hasattr(self.estimator, "coef_"):
            self.coef_ = self.estimator.coef_

        else:
            raise ValueError("Missing `feature_importances_` or `coef_`"
                             " attribute. ")

        return self
