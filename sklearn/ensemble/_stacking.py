from abc import ABCMeta
from copy import deepcopy

import numpy as np

from ..base import clone
from ..base import ClassifierMixin, RegressorMixin, TransformerMixin
from ..base import is_classifier, is_regressor
from ..base import MetaEstimatorMixin

from .base import _parallel_fit_estimator

from sklearn.externals.joblib import Parallel, delayed
from sklearn.externals.six import with_metaclass
from sklearn.externals.six import string_types

from sklearn.linear_model import LogisticRegression
from sklearn.linear_model import LinearRegression

from ..model_selection import cross_val_predict
from ..model_selection import check_cv

from ..utils import check_random_state
from ..utils.metaestimators import _BaseComposition
from ..utils.metaestimators import if_delegate_has_method
from ..utils.validation import has_fit_parameter
from ..utils.validation import check_is_fitted


class BaseStacking(with_metaclass(ABCMeta, _BaseComposition,
                                  MetaEstimatorMixin, TransformerMixin)):

    def __init__(self, base_estimators=None, stacking_estimator=None, cv=None,
                 method_base_estimators='auto', n_jobs=1, random_state=None,
                 verbose=0):
        self.base_estimators = base_estimators
        self.stacking_estimator = stacking_estimator
        self.cv = cv
        self.method_base_estimators = method_base_estimators
        self.n_jobs = n_jobs
        self.random_state = random_state
        self.verbose = verbose

    def _validate_stacking_estimator(self, default=None):
        if self.stacking_estimator is not None:
            self.stacking_estimator_ = clone(self.stacking_estimator)
        else:
            self.stacking_estimator_ = clone(default)

    @property
    def named_base_estimators(self):
        return dict(self.base_estimators)

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

        self._validate_meta_esimator()

        if self.base_estimators is None or len(self.base_estimators) == 0:
            raise AttributeError('Invalid `base_estimators` attribute, '
                                 '`base_estimators` should be a list of '
                                 '(string, estimator) tuples')

        if sample_weight is not None:
            for name, step in self.base_estimators:
                if not has_fit_parameter(step, 'sample_weight'):
                    raise ValueError('Underlying estimator \'%s\' does not'
                                     ' support sample weights.' % name)

        names, base_estimators_ = zip(*self.base_estimators)
        self._validate_names(names)

        if isinstance(self.method_base_estimators, string_types):
            if self.method_base_estimators != 'auto':
                raise AttributeError('When "method_base_estimators" is a '
                                     'string, it should be "auto". Got {} '
                                     'instead.'
                                     .format(self.method_base_estimators))
            method_base_estimators = ([self.method_base_estimators] *
                                      len(base_estimators_))
        else:
            if len(self.base_estimators) != len(self.method_base_estimators):
                raise AttributeError('When "method_base_estimators" is a '
                                     'list, it should be the same length as '
                                     'the list of base_estimators. Provided '
                                     '{} methods for {} base_estimators.'
                                     .format(len(self.method_base_estimators),
                                             len(self.base_estimators)))
            method_base_estimators = self.model_base_estimators

        self.method_base_estimators_ = [
            self._method_name(name, est, meth)
            for name, est, meth in zip(names, base_estimators_,
                                       method_base_estimators)]

        # Fit the base estimators on the whole training data. Those
        # base estimators will be used in transform, predict, and
        # predict_proba. There exposed publicly.
        self.base_estimators_ = Parallel(n_jobs=self.n_jobs)(
            delayed(_parallel_fit_estimator)(clone(est), X, y, sample_weight)
            for est in base_estimators_)

        # To train the meta-classifier using the most data as possible, we use
        # a cross-validation to predict the output of the stacked estimators.

        # To ensure that the data provided to each estimator are the same, we
        # need to set the random state of the cv if there is one and we need to
        # take a copy.
        random_state = check_random_state(self.random_state)
        cv = check_cv(self.cv)
        if hasattr(cv, 'random_state'):
            cv.random_state = random_state

        X_meta = Parallel(n_jobs=self.n_jobs)(
            delayed(cross_val_predict)(clone(est), X, y, cv=deepcopy(cv),
                                       method=meth, n_jobs=self.n_jobs,
                                       verbose=self.verbose)
            for est, meth in zip(base_estimators_,
                                 self.method_base_estimators_))
        X_meta = self._concatenate_predictions(X_meta)
        self.stacking_estimator_.fit(X_meta, y)

        return self

    def transform(self, X):
        check_is_fitted(self, 'base_estimators_')
        return self._concatenate_predictions([
            getattr(est, meth)(X)
            for est, meth in zip(self.base_estimators_,
                                 self.method_base_estimators_)])

    @if_delegate_has_method(delegate='stacking_estimator_')
    def predict(self, X):
        check_is_fitted(self, ['base_estimators_', 'stacking_estimator_'])
        return self.stacking_estimator_.predict(self.transform(X))

    @if_delegate_has_method(delegate='stacking_estimator_')
    def predict_proba(self, X):
        check_is_fitted(self, ['base_estimators_', 'stacking_estimator_'])
        return self.stacking_estimator_.predict_proba(self.transform(X))


class StackingClassifier(BaseStacking, ClassifierMixin):
    """Stacked of estimator using a classifier.

    Parameters
    ----------
    base_estimators : list of (string, estimator) tuples
        Base estimators which will be stacked together.

    stacking_estimator : estimator object
        A classifier which will be used to combine the base estimators.

    cv : int, cross-validation generator or an iterable, optional
        Determines the cross-validation splitting strategy. Possible inputs for
        cv are:

        * None, to use the default 3-fold cross validation,
        * integer, to specify the number of folds in a (Stratified) KFold,
        * An object to be used as a cross-validation generator,
        * An iterable yielding train, test splits.

        For integer/None inputs, if the estimator is a classifier and y is
        either binary or multiclass, StratifiedKFold is used. In all other
        cases, KFold is used.

        Refer :ref:`User Guide <cross_validation>` for the various
        cross-validation strategies that can be used here.

    method_base_estimators : list of string or 'auto', optional
        Methods called for each base estimator. It can be:

        * if a list of string in which each string is associated to the
          ``base_estimators``,
        * if ``auto``, it will try to invoke, for each estimator,
        ``predict_proba``, ``decision_function`` or ``predict`` in that order.

    n_jobs : int, optional (default=1)
        The number of jobs to ``fit`` the ``base_estimators`` in parallel. If
        -1, then the number of jobs is set to the number of cores.

    random_state : int, RandomState instance or None, optional (default=None)
        If int, random_state is the seed used by the random number generator;
        If RandomState instance, random_state is the random number generator;
        If None, the random number generator is the RandomState instance used
        by `np.random`. Used to set the ``cv``.

    Attributes
    ----------
    base_estimators_ : list of estimator object
        The base estimators fitted.

    stacking_classifier_ : estimator object
        The classifier to stacked the base estimators fitted.

    method_base_estimators_ : list of string
        The method used by each base estimator.

    Examples
    --------
    >>> from sklearn.datasets import load_iris
    >>> X, y = load_iris(return_X_y=True)
    >>> from sklearn.linear_model import LogisticRegression
    >>> from sklearn.svm import LinearSVC
    >>> from sklearn.ensemble import RandomForestClassifier
    >>> from sklearn.ensemble import StackingClassifier
    >>> base_estimators = [('lr', LogisticRegression()), ('svr', LinearSVC())]
    >>> clf = StackingClassifier(base_estimators=base_estimators,
    ...                          stacking_estimator=RandomForestClassifier())
    >>> from sklearn.model_selection import train_test_split
    >>> X_train, X_test, y_train, y_test = train_test_split(X, y)
    >>> clf.fit(X_train, y_train).score(X_test, y_test) # doctest: +ELLIPSIS
    0...


    """
    def __init__(self, base_estimators=None, stacking_estimator=None, cv=None,
                 method_base_estimators='auto', n_jobs=1, random_state=None,
                 verbose=0):
        super(StackingClassifier, self).__init__(
            base_estimators=base_estimators,
            stacking_estimator=stacking_estimator,
            cv=cv,
            method_base_estimators=method_base_estimators,
            n_jobs=n_jobs,
            random_state=random_state,
            verbose=verbose)

    def _validate_meta_esimator(self):
        super(StackingClassifier, self)._validate_stacking_estimator(
            default=LogisticRegression(random_state=self.random_state))
        if not is_classifier(self.stacking_estimator_):
            raise AttributeError('`stacking_estimator` attribute should be a '
                                 'classifier.')


class StackingRegressor(BaseStacking, RegressorMixin):
    """Stacked of estimator using a regressor.

    Parameters
    ----------
    base_estimators : list of (string, estimator) tuples
        Base estimators which will be stacked together.

    stacking_estimator : estimator object
        A regressor which will be used to combine the base estimators.

    cv : int, cross-validation generator or an iterable, optional
        Determines the cross-validation splitting strategy. Possible inputs for
        cv are:

        * None, to use the default 3-fold cross validation,
        * integer, to specify the number of folds in a (Stratified) KFold,
        * An object to be used as a cross-validation generator,
        * An iterable yielding train, test splits.

        For integer/None inputs, if the estimator is a classifier and y is
        either binary or multiclass, StratifiedKFold is used. In all other
        cases, KFold is used.

        Refer :ref:`User Guide <cross_validation>` for the various
        cross-validation strategies that can be used here.

    method_base_estimators : list of string or 'auto', optional
        Methods called for each base estimator. It can be:

        * if a list of string in which each string is associated to the
          ``base_estimators``,
        * if ``auto``, it will try to invoke, for each estimator,
        ``predict_proba``, ``decision_function`` or ``predict`` in that order.

    n_jobs : int, optional (default=1)
        The number of jobs to ``fit`` the ``base_estimators`` in parallel. If
        -1, then the number of jobs is set to the number of cores.

    random_state : int, RandomState instance or None, optional (default=None)
        If int, random_state is the seed used by the random number generator;
        If RandomState instance, random_state is the random number generator;
        If None, the random number generator is the RandomState instance used
        by `np.random`. Used to set the ``cv``.

    Attributes
    ----------
    base_estimators_ : list of estimator object
        The base estimators fitted.

    stacking_estimator_ : estimator object
        The regressor to stacked the base estimators fitted.

    method_base_estimators_ : list of string
        The method used by each base estimator.

    Examples
    --------
    >>> from sklearn.datasets import load_diabetes
    >>> X, y = load_diabetes(return_X_y=True)
    >>> from sklearn.linear_model import LinearRegression
    >>> from sklearn.svm import LinearSVR
    >>> from sklearn.ensemble import RandomForestRegressor
    >>> from sklearn.ensemble import StackingRegressor
    >>> base_estimators = [('lr', LinearRegression()), ('svr', LinearSVR())]
    >>> clf = StackingRegressor(base_estimators=base_estimators,
    ...                         stacking_estimator=RandomForestRegressor())
    >>> from sklearn.model_selection import train_test_split
    >>> X_train, X_test, y_train, y_test = train_test_split(X, y)
    >>> clf.fit(X_train, y_train).score(X_test, y_test) # doctest: +ELLIPSIS
    0...

    """
    def __init__(self, base_estimators=None, stacking_estimator=None, cv=None,
                 method_base_estimators='auto', n_jobs=1, random_state=None,
                 verbose=0):
        super(StackingRegressor, self).__init__(
            base_estimators=base_estimators,
            stacking_estimator=stacking_estimator,
            cv=cv,
            method_base_estimators=method_base_estimators,
            n_jobs=n_jobs,
            random_state=random_state,
            verbose=verbose)

    def _validate_meta_esimator(self):
        super(StackingRegressor, self)._validate_stacking_estimator(
            default=LinearRegression())
        if not is_regressor(self.stacking_estimator_):
            raise AttributeError('`stacking_estimator` attribute should be a '
                                 'regressor.')
