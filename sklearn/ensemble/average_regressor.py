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
from ..base import clone
from ..externals.joblib import Parallel, delayed
from ..utils import check_array
from ..utils.metaestimators import _BaseComposition
from ..utils.validation import check_is_fitted

from .base import BaseEnsemble, _partition_estimators


__all__ = ["AverageRegressor"]


def _parallel_fit_estimator(estimator, X, y, sample_weight=None):
    """Private function used to fit an estimator within a job."""
    if sample_weight is not None:
        estimator.fit(X, y, sample_weight=sample_weight)
    else:
        estimator.fit(X, y)
    return estimator


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

    Read more in the :ref:`User Guide <average_regressor>`.

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

    def __init__(self, estimators, weights=None, verbose=0, n_jobs=1):
        self.estimators = estimators
        self.weights = weights
        self.n_jobs = n_jobs
        self.verbose = verbose

    @property
    def named_estimators(self):
        return Bunch(**dict(self.estimators))

    @property
    def _weights_not_none(self):
        """Get the weights of not `None` estimators"""
        if self.weights is None:
            return None
        return [w for est, w in zip(self.estimators,
                                    self.weights) if est[1] is not None]

    def fit(self, X, y, sample_weight=None):
        """ Fit the estimators.

        Parameters
        ----------
        X : {array-like, sparse matrix}, shape = [n_samples, n_features]
            Training vectors, where n_samples is the number of samples and
            n_features is the number of features.

        y : array-like, shape = [n_samples]
            Target values.

        sample_weight : array-like, shape = [n_samples] or None
            Sample weights. If None, then samples are equally weighted.
            Note that this is supported only if all underlying estimators
            support sample weights.

        Returns
        -------
        self : object
        """

        # validate that the estimaors is not an empty list
        if self.estimators is None or len(self.estimators) == 0:
            raise AttributeError('Invalid `estimators` attribute, `estimators`'
                                 ' should be a list of (string, estimator)'
                                 ' tuples')

        # validate that at least one estimator is not None
        isnone = [clf is not None for _, clf in self.estimators]
        if not any(isnone):
            raise ValueError('All estimators are None. At least one is '
                             'required to be a regressor!')

        # check that len(weights) == len(estimators)
        if (self.weights is not None and
                len(self.weights) != len(self.estimators)):
            raise ValueError('Number of classifiers and weights must be equal'
                             '; got %d weights, %d estimators'
                             % (len(self.weights), len(self.estimators)))

        names, clfs = zip(*self.estimators)
        self._validate_names(names)

        # fit estimators in parallel if n_jobs > 1
        self.estimators_ = Parallel(n_jobs=self.n_jobs)(
                delayed(_parallel_fit_estimator)(clone(clf), X, y,
                                                 sample_weight=sample_weight)
                for clf in clfs if clf is not None)

        return self

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
        n_estimators = len(self.estimators_) 
        n_jobs, n_estimators, starts = _partition_estimators(n_estimators,
                                                             self.n_jobs)

        all_y_hat = Parallel(n_jobs=n_jobs, verbose=self.verbose)(
            delayed(_parallel_predict_regression)(
                self.estimators_[starts[i]:starts[i + 1]],
                self.estimators_features_[starts[i]:starts[i + 1]],
                X)
            for i in range(n_jobs))

        # Reduce
        y_hat = np.average(all_y_hat, axis=0,
                           weights=self._weights_not_none)
        
        return y_hat
