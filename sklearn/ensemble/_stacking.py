from copy import deepcopy

import numpy as np

from ..base import clone
from ..base import ClassifierMixin, RegressorMixin, TransformerMixin
from ..base import is_classifier, is_regressor
from ..base import MetaEstimatorMixin

from .base import _parallel_fit_estimator

from sklearn.externals.joblib import Parallel, delayed
from sklearn.externals.six import string_types

from sklearn.linear_model import LogisticRegression

from ..model_selection import cross_val_predict
from ..model_selection import check_cv

from ..utils import check_random_state
from ..utils.metaestimators import _BaseComposition
from ..utils.metaestimators import if_delegate_has_method
from ..utils.validation import has_fit_parameter
from ..utils.validation import check_is_fitted


class StackingClassifier(_BaseComposition, MetaEstimatorMixin, ClassifierMixin,
                         TransformerMixin):

    def __init__(self, estimators=None, classifier=None, cv=3,
                 method_estimators='auto', n_jobs=1, random_state=None,
                 verbose=0):
        self.estimators = estimators
        self.classifier = classifier
        self.cv = cv
        self.method_estimators = method_estimators
        self.n_jobs = n_jobs
        self.random_state = random_state
        self.verbose = verbose

    @property
    def named_estimators(self):
        return dict(self.estimators)

    @staticmethod
    def _method_name(name, estimator, method):
        if method == 'auto':
            if getattr(estimator, 'predict_proba', None):
                return 'predict_proba'
            elif getattr(estimator, 'decision_function', None):
                return 'decision_function'
            else:
                return 'predict'
        else:
            if not hasattr(estimator, method):
                raise ValueError('Underlying estimator {} does not implement '
                                 'the method {}.'.format(name, method))
            return method

    @staticmethod
    def _concatenate_predictions(y_pred):
        return np.concatenate([pred.reshape(-1, 1) if pred.ndim == 1
                               else pred for pred in y_pred], axis=1)

    def fit(self, X, y, sample_weight=None):

        random_state = check_random_state(self.random_state)

        if self.classifier is None:
            self.classifier_ = LogisticRegression(random_state=random_state)
        else:
            if not is_classifier(self.classifier):
                raise AttributeError('`classifier` attribute should be a '
                                     'classifier.')
            self.classifier_ = clone(self.classifier)

        if self.estimators is None or len(self.estimators) == 0:
            raise AttributeError('Invalid `estimators` attribute, `estimators`'
                                 ' should be a list of (string, estimator)'
                                 ' tuples')

        if sample_weight is not None:
            for name, step in self.estimators:
                if not has_fit_parameter(step, 'sample_weight'):
                    raise ValueError('Underlying estimator \'%s\' does not'
                                     ' support sample weights.' % name)

        names, estimators_ = zip(*self.estimators)
        self._validate_names(names)

        if isinstance(self.method_estimators, string_types):
            if self.method_estimators != 'auto':
                raise AttributeError('When "method" is a string, it should be '
                                     '"auto". Got {} instead.'
                                     .format(self.method))
            method_estimators = [self.method_estimators] * len(estimators_)
        else:
            if len(self.estimators) != len(self.method_estimators):
                raise AttributeError('When "method" is a list, it should be '
                                     'the same length as the list of '
                                     'estimators. Provided {} methods for {} '
                                     'estimators.'
                                     .format(len(self.method_estimators),
                                             len(self.estimators)))
            method_estimators = self.model_estimators

        self.method_estimators_ = [
            self._method_name(name, est, meth)
            for name, est, meth in zip(names, estimators_, method_estimators)]

        # Fit the estimators on the whole training data. Those estimators will
        # be used in transform, predict, and predict_proba. There exposed
        # publicly.
        self.estimators_ = Parallel(n_jobs=self.n_jobs)(
            delayed(_parallel_fit_estimator)(clone(est), X, y, sample_weight)
            for est in estimators_)

        # To train the meta-classifier using the most data as possible, we use
        # a cross-validation to predict the output of the stacked estimators.

        # To ensure that the data provided to each estimator are the same, we
        # need to set the random state of the cv if there is one and we need to
        # take a copy.
        cv = check_cv(self.cv)
        if hasattr(cv, 'random_state'):
            cv.random_state = random_state

        X_meta = Parallel(n_jobs=self.n_jobs)(
            delayed(cross_val_predict)(clone(est), X, y, cv=deepcopy(cv),
                                       method=meth, n_jobs=self.n_jobs,
                                       verbose=self.verbose)
            for est, meth in zip(estimators_, self.method_estimators_))
        X_meta = self._concatenate_predictions(X_meta)
        self.classifier_.fit(X_meta, y)

        return self

    def transform(self, X):
        check_is_fitted(self, 'estimators_')
        return self._concatenate_predictions([
            getattr(est, meth)(X)
            for est, meth in zip(self.estimators_, self.method_estimators_)])

    @if_delegate_has_method(delegate='classifier_')
    def predict(self, X):
        check_is_fitted(self, ['estimators_', 'classifier_'])
        return self.classifier_.predict(self.transform(X))

    @if_delegate_has_method(delegate='classifier_')
    def predict_proba(self, X):
        check_is_fitted(self, ['estimators_', 'classifier_'])
        return self.classifier_.predict_proba(self.transform(X))
