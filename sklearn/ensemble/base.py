"""
Base class for ensemble-based estimators.
"""

# Authors: Gilles Louppe
# License: BSD 3 clause

import numpy as np

from ..base import clone
from ..base import BaseEstimator
from ..base import MetaEstimatorMixin
from ..externals.joblib import cpu_count


class BaseEnsemble(BaseEstimator, MetaEstimatorMixin):
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

    Attributes
    ----------
    base_estimator_ : list of estimators
        The base estimator from which the ensemble is grown.

    estimators_ : list of estimators
        The collection of fitted base estimators.
    """

    def __init__(self, base_estimator, n_estimators=10,
                 estimator_params=tuple()):
        # Set parameters
        self.base_estimator = base_estimator
        self.n_estimators = n_estimators
        self.estimator_params = estimator_params

        # Don't instantiate estimators now! Parameters of base_estimator might
        # still change. Eg., when grid-searching with the nested object syntax.
        # This needs to be filled by the derived classes.
        self.estimators_ = []

    def _validate_estimator(self, default=None):
        """Check the estimator and the n_estimator attribute, set the
        `base_estimator_` attribute."""
        if self.n_estimators <= 0:
            raise ValueError("n_estimators must be greater than zero, "
                             "got {0}.".format(self.n_estimators))

        if self.base_estimator is not None:
            self.base_estimator_ = self.base_estimator
        else:
            self.base_estimator_ = default

        if self.base_estimator_ is None:
            raise ValueError("base_estimator cannot be None")

    def _make_estimator(self, append=True):
        """Make and configure a copy of the `base_estimator_` attribute.

        Warning: This method should be used to properly instantiate new
        sub-estimators.
        """
        estimator = clone(self.base_estimator_)
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

    def __iter__(self):
        """Returns iterator over estimators in the ensemble."""
        return iter(self.estimators_)


def _partition_estimators(n_estimators, n_jobs):
    """Private function used to partition estimators between jobs."""
    # Compute the number of jobs
    if n_jobs == -1:
        n_jobs = min(cpu_count(), n_estimators)

    else:
        n_jobs = min(n_jobs, n_estimators)

    # Partition estimators between jobs
    n_estimators_per_job = (n_estimators // n_jobs) * np.ones(n_jobs,
                                                              dtype=np.int)
    n_estimators_per_job[:n_estimators % n_jobs] += 1
    starts = np.cumsum(n_estimators_per_job)

    return n_jobs, n_estimators_per_job.tolist(), [0] + starts.tolist()
