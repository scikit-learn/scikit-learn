"""
This module implements multioutput regression and classification.

The estimators provided in this module are meta-estimators: they require
a base estimator to be provided in their constructor. The meta-estimator
extends single output estimators to multioutput estimators.
"""

# Author: Tim Head <betatim@gmail.com>
# Author: Hugo Bowne-Anderson <hugobowne@gmail.com>
# Author: Chris Rivera <chris.richard.rivera@gmail.com>
# Author: Michael Williamson
# Author: James Ashton Nichols <james.ashton.nichols@gmail.com>
#
# License: BSD 3 clause

import numpy as np

from abc import ABCMeta
from .base import BaseEstimator, clone
from .base import RegressorMixin, ClassifierMixin
from .utils import check_array, check_X_y
from .utils.fixes import parallel_helper
from .utils.validation import check_is_fitted, has_fit_parameter
from .utils.metaestimators import if_delegate_has_method
from .externals.joblib import Parallel, delayed
from .externals import six

__all__ = ["MultiOutputRegressor", "MultiOutputClassifier"]


def _fit_estimator(estimator, X, y, sample_weight=None):
    estimator = clone(estimator)
    if sample_weight is not None:
        estimator.fit(X, y, sample_weight=sample_weight)
    else:
        estimator.fit(X, y)
    return estimator


def _partial_fit_estimator(estimator, X, y, classes=None, sample_weight=None,
                           first_time=True):
    if first_time:
        estimator = clone(estimator)

    if sample_weight is not None:
        if classes is not None:
            estimator.partial_fit(X, y, classes=classes,
                                  sample_weight=sample_weight)
        else:
            estimator.partial_fit(X, y, sample_weight=sample_weight)
    else:
        if classes is not None:
            estimator.partial_fit(X, y, classes=classes)
        else:
            estimator.partial_fit(X, y)
    return estimator


class MultiOutputEstimator(six.with_metaclass(ABCMeta, BaseEstimator)):

    def __init__(self, estimator, n_jobs=1):
        self.estimator = estimator
        self.n_jobs = n_jobs

    @if_delegate_has_method('estimator')
    def partial_fit(self, X, y, classes=None, sample_weight=None):
        """Incrementally fit the model to data.
        Fit a separate model for each output variable.

        Parameters
        ----------
        X : (sparse) array-like, shape (n_samples, n_features)
            Data.

        y : (sparse) array-like, shape (n_samples, n_outputs)
            Multi-output targets.

        classes : list of numpy arrays, shape (n_outputs)
            Each array is unique classes for one output in str/int
            Can be obtained by via
            ``[np.unique(y[:, i]) for i in range(y.shape[1])]``, where y is the
            target matrix of the entire dataset.
            This argument is required for the first call to partial_fit
            and can be omitted in the subsequent calls.
            Note that y doesn't need to contain all labels in `classes`.

        sample_weight : array-like, shape = (n_samples) or None
            Sample weights. If None, then samples are equally weighted.
            Only supported if the underlying regressor supports sample
            weights.

        Returns
        -------
        self : object
            Returns self.
        """
        X, y = check_X_y(X, y,
                         multi_output=True,
                         accept_sparse=True)

        if y.ndim == 1:
            raise ValueError("y must have at least two dimensions for "
                             "multi-output regression but has only one.")

        if (sample_weight is not None and
                not has_fit_parameter(self.estimator, 'sample_weight')):
            raise ValueError("Underlying estimator does not support"
                             " sample weights.")

        first_time = not hasattr(self, 'estimators_')

        self.estimators_ = Parallel(n_jobs=self.n_jobs)(
            delayed(_partial_fit_estimator)(
                self.estimators_[i] if not first_time else self.estimator,
                X, y[:, i],
                classes[i] if classes is not None else None,
                sample_weight, first_time) for i in range(y.shape[1]))
        return self

    def fit(self, X, y, sample_weight=None):
        """ Fit the model to data.
        Fit a separate model for each output variable.

        Parameters
        ----------
        X : (sparse) array-like, shape (n_samples, n_features)
            Data.

        y : (sparse) array-like, shape (n_samples, n_outputs)
            Multi-output targets. An indicator matrix turns on multilabel
            estimation.

        sample_weight : array-like, shape = (n_samples) or None
            Sample weights. If None, then samples are equally weighted.
            Only supported if the underlying regressor supports sample
            weights.

        Returns
        -------
        self : object
            Returns self.
        """

        if not hasattr(self.estimator, "fit"):
            raise ValueError("The base estimator should implement a fit method")

        X, y = check_X_y(X, y,
                         multi_output=True,
                         accept_sparse=True)

        if y.ndim == 1:
            raise ValueError("y must have at least two dimensions for "
                             "multi-output regression but has only one.")

        if (sample_weight is not None and
                not has_fit_parameter(self.estimator, 'sample_weight')):
            raise ValueError("Underlying estimator does not support"
                             " sample weights.")

        self.estimators_ = Parallel(n_jobs=self.n_jobs)(
            delayed(_fit_estimator)(
                self.estimator, X, y[:, i], sample_weight)
            for i in range(y.shape[1]))
        return self

    def predict(self, X):
        """Predict multi-output variable using a model
         trained for each target variable.

        Parameters
        ----------
        X : (sparse) array-like, shape (n_samples, n_features)
            Data.

        Returns
        -------
        y : (sparse) array-like, shape (n_samples, n_outputs)
            Multi-output targets predicted across multiple predictors.
            Note: Separate models are generated for each predictor.
        """
        check_is_fitted(self, 'estimators_')
        if not hasattr(self.estimator, "predict"):
            raise ValueError("The base estimator should implement a predict method")

        X = check_array(X, accept_sparse=True)

        y = Parallel(n_jobs=self.n_jobs)(
            delayed(parallel_helper)(e, 'predict', X)
            for e in self.estimators_)

        return np.asarray(y).T


class MultiOutputRegressor(MultiOutputEstimator, RegressorMixin):
    """Multi target regression

    This strategy consists of fitting one regressor per target. This is a
    simple strategy for extending regressors that do not natively support
    multi-target regression.

    Parameters
    ----------
    estimator : estimator object
        An estimator object implementing `fit` and `predict`.

    n_jobs : int, optional, default=1
        The number of jobs to run in parallel for `fit`. If -1,
        then the number of jobs is set to the number of cores.
        When individual estimators are fast to train or predict
        using `n_jobs>1` can result in slower performance due
        to the overhead of spawning processes.
    """

    def __init__(self, estimator, n_jobs=1):
        super(MultiOutputRegressor, self).__init__(estimator, n_jobs)

    def partial_fit(self, X, y, sample_weight=None):
        """Incrementally fit the model to data.
        Fit a separate model for each output variable.

        Parameters
        ----------
        X : (sparse) array-like, shape (n_samples, n_features)
            Data.

        y : (sparse) array-like, shape (n_samples, n_outputs)
            Multi-output targets.

        sample_weight : array-like, shape = (n_samples) or None
            Sample weights. If None, then samples are equally weighted.
            Only supported if the underlying regressor supports sample
            weights.

        Returns
        -------
        self : object
            Returns self.
        """
        super(MultiOutputRegressor, self).partial_fit(
            X, y, sample_weight=sample_weight)

    def score(self, X, y, sample_weight=None):
        """Returns the coefficient of determination R^2 of the prediction.

        The coefficient R^2 is defined as (1 - u/v), where u is the regression
        sum of squares ((y_true - y_pred) ** 2).sum() and v is the residual
        sum of squares ((y_true - y_true.mean()) ** 2).sum().
        Best possible score is 1.0 and it can be negative (because the
        model can be arbitrarily worse). A constant model that always
        predicts the expected value of y, disregarding the input features,
        would get a R^2 score of 0.0.

        Notes
        -----
        R^2 is calculated by weighting all the targets equally using
        `multioutput='uniform_average'`.

        Parameters
        ----------
        X : array-like, shape (n_samples, n_features)
            Test samples.

        y : array-like, shape (n_samples) or (n_samples, n_outputs)
            True values for X.

        sample_weight : array-like, shape [n_samples], optional
            Sample weights.

        Returns
        -------
        score : float
            R^2 of self.predict(X) wrt. y.
        """
        # XXX remove in 0.19 when r2_score default for multioutput changes
        from .metrics import r2_score
        return r2_score(y, self.predict(X), sample_weight=sample_weight,
                        multioutput='uniform_average')


class MultiOutputClassifier(MultiOutputEstimator, ClassifierMixin):
    """Multi target classification

    This strategy consists of fitting one classifier per target. This is a
    simple strategy for extending classifiers that do not natively support
    multi-target classification

    Parameters
    ----------
    estimator : estimator object
        An estimator object implementing `fit`, `score` and `predict_proba`.

    n_jobs : int, optional, default=1
        The number of jobs to use for the computation. If -1 all CPUs are used.
        If 1 is given, no parallel computing code is used at all, which is
        useful for debugging. For n_jobs below -1, (n_cpus + 1 + n_jobs) are
        used. Thus for n_jobs = -2, all CPUs but one are used.
        The number of jobs to use for the computation.
        It does each target variable in y in parallel.

    Attributes
    ----------
    estimators_ : list of `n_output` estimators
        Estimators used for predictions.
    """

    def __init__(self, estimator, n_jobs=1):
        super(MultiOutputClassifier, self).__init__(estimator, n_jobs)

    def predict_proba(self, X):
        """Probability estimates.
        Returns prediction probabilites for each class of each output.

        Parameters
        ----------
        X : array-like, shape (n_samples, n_features)
            Data

        Returns
        -------
        p : array of shape = [n_samples, n_classes], or a list of n_outputs \
            such arrays if n_outputs > 1.
            The class probabilities of the input samples. The order of the
            classes corresponds to that in the attribute `classes_`.
        """
        check_is_fitted(self, 'estimators_')
        if not hasattr(self.estimator, "predict_proba"):
            raise ValueError("The base estimator should implement"
                             "predict_proba method")

        results = [estimator.predict_proba(X) for estimator in
                   self.estimators_]
        return results

    def score(self, X, y):
        """"Returns the mean accuracy on the given test data and labels.

        Parameters
        ----------
        X : array-like, shape [n_samples, n_features]
            Test samples

        y : array-like, shape [n_samples, n_outputs]
            True values for X

        Returns
        -------
        scores : float
            accuracy_score of self.predict(X) versus y
        """
        check_is_fitted(self, 'estimators_')
        n_outputs_ = len(self.estimators_)
        if y.ndim == 1:
            raise ValueError("y must have at least two dimensions for "
                             "multi target classification but has only one")
        if y.shape[1] != n_outputs_:
            raise ValueError("The number of outputs of Y for fit {0} and"
                             " score {1} should be same".
                             format(n_outputs_, y.shape[1]))
        y_pred = self.predict(X)
        return np.mean(np.all(y == y_pred, axis=1))
