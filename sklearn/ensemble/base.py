"""
Base class for ensemble-based estimators.
"""

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

    **estimator_params : key-words parameters
        The parameters to set for the underlying estimator. If none are given,
        default parameters are used.
    """
    def __init__(self, base_estimator, n_estimators, **estimator_params):
        # Check parameters
        if not isinstance(base_estimator, BaseEstimator):
            raise TypeError("estimator must be a subclass of BaseEstimator")
        if n_estimators <= 0:
            raise ValueError("n_estimators must be greater than zero.")

        # Make clones
        self.base_estimator = base_estimator
        self.base_estimator.set_params(**estimator_params)
        self.n_estimators = n_estimators
        self.estimators = [clone(base_estimator) for i in xrange(n_estimators)]

    def __len__(self):
        """Returns the number of estimators in the ensemble."""
        return self.n_estimators

    def __getitem__(self, index):
        """Returns the index'th estimator in the ensemble."""
        return self.estimators[index]
