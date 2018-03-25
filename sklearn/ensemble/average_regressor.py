"""
Average regressor.

This module contains an average regressor for
regression estimators.

"""

# Author: Mohamed Ali Jamaoui <m.ali.jamaoui@gmail.com>
# License: BSD 3 clause

import numpy as np

from ..base import RegressorMixin
from ..base import TransformerMixin
from ..base import clone
from ..externals.joblib import Parallel, delayed
from ..utils import check_array
from ..utils.metaestimators import _BaseComposition
from ..utils.validation import check_is_fitted, has_fit_parameter
from ..utils import Bunch


__all__ = ["AverageRegressor"]


def _parallel_fit_estimator(estimator, X, y, sample_weight=None):
    """Private function used to fit an estimator within a job."""
    if sample_weight is not None:
        estimator.fit(X, y, sample_weight=sample_weight)
    else:
        estimator.fit(X, y)
    return estimator


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

    verbose : int, optional (default=0)
        Controls the verbosity of the building process.

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
            raise ValueError('Number of regressors and weights must be equal'
                             '; got %d weights, %d estimators'
                             % (len(self.weights), len(self.estimators)))

        if sample_weight is not None:
            for name, step in self.estimators:
                if not has_fit_parameter(step, 'sample_weight'):
                    raise ValueError('Underlying estimator \'%s\' does not'
                                     ' support sample weights.' % name)

        names, clfs = zip(*self.estimators)
        self._validate_names(names)

        # fit estimators in parallel if n_jobs > 1
        self.estimators_ = Parallel(n_jobs=self.n_jobs)(
                delayed(_parallel_fit_estimator)(clone(clf), X, y,
                                                 sample_weight=sample_weight)
                for clf in clfs if clf is not None)

        self.named_estimators_ = Bunch(**dict())
        for k, e in zip(self.estimators, self.estimators_):
            self.named_estimators_[k[0]] = e

        return self

    def _collect_predictions(self, X):
        """Collect results from reg.predict calls. """
        return np.asarray([reg.predict(X) for reg in self.estimators_])

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

        y_hat = np.average(self._collect_predictions(X), axis=0,
                           weights=self._weights_not_none)
        return y_hat

    def set_params(self, **params):
        """ Setting the parameters for the average regressor

        Valid parameter keys can be listed with get_params().

        Parameters
        ----------
        params : keyword arguments
            Specific parameters using e.g. set_params(parameter_name=new_value)
            In addition, to setting the parameters of the ``AverageRegressor``,
            the individual regressors of the ``AverageRegressor`` can also be
            set or replaced by setting them to None.

        Examples
        --------
        >># TODO
        """
        super(AverageRegressor, self)._set_params('estimators', **params)
        return self

    def get_params(self, deep=True):
        """ Get the parameters of the AverageRegressor

        Parameters
        ----------
        deep: bool
            Setting it to True gets the various regressors and their parameters
        """
        return super(AverageRegressor,
                     self)._get_params('estimators', deep=deep)
