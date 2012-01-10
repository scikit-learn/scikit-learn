"""
Base class for ensemble-based estimators.
"""

# Authors: Gilles Louppe
# License: BSD 3

from ..base import clone
from ..base import BaseEstimator


class BaseEnsemble(BaseEstimator):
    """Base class for all ensemble classes.

    Warning: This class should not be used directly. Use derived classes
    instead.

    Parameters
    ----------
    base_estimator : object, optional (default=None)
        The base estimator from which the ensemble is built.

    n_estimators : integer
        The number of estimators in the ensemble.

    estimator_params : list of strings
        The list of attributes to use as parameters when instantiating a
        new base estimator. If none are given, default parameters are used.
    """
    def __init__(self, base_estimator, n_estimators, estimator_params=[]):
        # Check parameters
        if not isinstance(base_estimator, BaseEstimator):
            raise TypeError("estimator must be a subclass of BaseEstimator")
        if n_estimators <= 0:
            raise ValueError("n_estimators must be greater than zero.")

        # Set parameters
        self.base_estimator = base_estimator
        self.n_estimators = n_estimators
        self.estimator_params = estimator_params

        # Don't instantiate estimators now! Parameters of base_estimator might
        # still change. Eg., when grid-searching with the nested object syntax.
        # This needs to be filled by the derived classes.
        self.estimators_ = []

    def _make_estimator(self, append=True):
        """Makes, configures and returns a copy of the base estimator.

        Warning: This method should be used to properly instantiate new
        sub-estimators.
        """
        estimator = clone(self.base_estimator)
        estimator.set_params(**dict((p, getattr(self, p))
                                    for p in self.estimator_params))

        if append:
            self.estimators_.append(estimator)

        return estimator

    def __len__(self):
        """Returns the number of estimators in the ensemble."""
        return len(self.estimators_)

    def __getitem__(self, index):
        """Returns the index'th estimator in the ensemble."""
        return self.estimators_[index]
