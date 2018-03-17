"""
Average regressor.

This module contains an average regressor for
regression estimators.

"""

# Author: Mohamed Ali Jamaoui <m.ali.jamaoui@gmail.com>
#
# License: BSD 3 clause

from ..base import RegressorMixin
from ..base import TransformerMixin
from ..externals.joblib import Parallel, delayed
from ..utils import check_array
from ..utils.metaestimators import _BaseComposition
from ..utils.validation import check_is_fitted

from .base import BaseEnsemble, _partition_estimators


__all__ = ["AverageRegressor"]


def _parallel_predict_regression(estimators, estimators_features, X):
    """Private function used to compute predictions within a job."""
    return sum(estimator.predict(X[:, features])
               for estimator, features in zip(estimators,
                                              estimators_features))


class AverageRegressor(_BaseComposition, RegressorMixin, TransformerMixin):
    """
    An average regressor is an ensemble meta-estimator that fits base
    regressors each on the whole dataset. it, then, averages the individual
    predictions to form a final prediction.

        .. versionadded:: 0.20

    Read more in the :ref:`User Guide <voting_classifier>`.

    Parameters
    ----------
    estimators : list of (string, estimator) tuples
        Invoking the ``fit`` method on the ``AverageRegressor`` will fit clones
        of those original estimators that will be stored in the class attribute
        ``self.estimators_``. An estimator can be set to `None` using
        ``set_params``.

    weights : array-like, shape = [n_regressors], optional (default=`None`)
        Sequence of weights (`float` or `int`) to weight the prediction of 
        a given estimator. 

    n_jobs : int, optional (default=1)
        The number of jobs to run in parallel for ``fit``.
        If -1, then the number of jobs is set to the number of cores.

    Attributes
    ----------
    estimators_ : list of regressors
        The collection of fitted sub-estimators as defined in ``estimators``
        that are not `None`.

    named_estimators_ : Bunch object, a dictionary with attribute access
        Attribute to access any fitted sub-estimators by name.

        .. versionadded:: 0.20

    Examples
    --------
    >>> import numpy as np
    >>> # TODO 
    """

    def __init__(self, estimators, weights=None, n_jobs=1):
        self.estimators = estimators
        self.weights = weights
        self.n_jobs = n_jobs

    def _set_oob_score(self, X, y):
        pass

    def fit(self):
        pass

    def predict(self, X):
        """Predict regression target for X.

        The predicted regression target of an input sample is computed as the
        mean predicted regression targets of the estimators in the ensemble.

        Parameters
        ----------
        X : {array-like, sparse matrix} of shape = [n_samples, n_features]
            The training input samples. Sparse matrices are accepted only if
            they are supported by the base estimator.

        Returns
        -------
        y : array of shape = [n_samples]
            The predicted values.
        """
        check_is_fitted(self, "estimators_")
        # Check data
        X = check_array(X, accept_sparse=['csr', 'csc'])

        # Parallel loop
        n_jobs, n_estimators, starts = _partition_estimators(self.n_estimators,
                                                             self.n_jobs)

        all_y_hat = Parallel(n_jobs=n_jobs, verbose=self.verbose)(
            delayed(_parallel_predict_regression)(
                self.estimators_[starts[i]:starts[i + 1]],
                self.estimators_features_[starts[i]:starts[i + 1]],
                X)
            for i in range(n_jobs))

        # Reduce
        y_hat = sum(all_y_hat) / self.n_estimators

        return y_hat

    def transform():
        pass

    def predict_proba():
        pass
