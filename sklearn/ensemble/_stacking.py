"""Stacking classifier and regressor."""

# Authors: Guillaume Lemaitre <g.lemaitre58@gmail.com>
# License: BSD 3 clause

from abc import ABCMeta, abstractmethod
from copy import deepcopy

import numpy as np
from joblib import Parallel, delayed

from ..base import clone
from ..base import ClassifierMixin, RegressorMixin, TransformerMixin
from ..base import is_classifier, is_regressor
from ..base import MetaEstimatorMixin

from .base import _parallel_fit_estimator

from ..linear_model import LogisticRegression
from ..linear_model import LinearRegression

from ..model_selection import cross_val_predict
from ..model_selection import check_cv

from ..utils import check_random_state
from ..utils import Bunch
from ..utils.metaestimators import _BaseComposition
from ..utils.metaestimators import if_delegate_has_method
from ..utils.validation import check_is_fitted


class _BaseStacking(_BaseComposition, MetaEstimatorMixin, TransformerMixin,
                    metaclass=ABCMeta):
    """Base class for stacking method.

    Warning: This class should not be used directly. Use derived classes
    instead.
    """
    _required_parameters = ['estimators']

    @abstractmethod
    def __init__(self, estimators, final_estimator=None, cv=None,
                 predict_method='auto', n_jobs=None, random_state=None,
                 verbose=0):
        self.estimators = estimators
        self.final_estimator = final_estimator
        self.cv = cv
        self.predict_method = predict_method
        self.n_jobs = n_jobs
        self.random_state = random_state
        self.verbose = verbose

    def _validate_final_estimator(self, default=None):
        if self.final_estimator is not None:
            self.final_estimator_ = clone(self.final_estimator)
        else:
            self.final_estimator_ = clone(default)

    @staticmethod
    def _method_name(name, estimator, method):
        if estimator in (None, 'drop'):
            return None
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

    def set_params(self, **params):
        """Setting the parameters for the stacking estimator.

        Valid parameter keys can be listed with get_params().

        Parameters
        ----------
        params : keyword arguments
            Specific parameters using e.g. set_params(parameter_name=new_value)
            In addition, to setting the parameters of the stacking estimator,
            the individual estimator of the stacking estimator can also be
            set or replaced by setting them to None.

        Examples
        --------
        # In this example, the RandomForestClassifier is removed
        clf1 = LogisticRegression()
        clf2 = RandomForestClassifier()
        eclf = StackingClassifier(estimators=[('lr', clf1), ('rf', clf2)]
        eclf.set_params(rf=None)

        """
        super()._set_params('estimators', **params)
        return self

    def get_params(self, deep=True):
        """ Get the parameters of the stacking estimator.

        Parameters
        ----------
        deep: bool
            Setting it to True gets the various classifiers and the parameters
            of the classifiers as well.
        """
        return super()._get_params('estimators', deep=deep)

    def _concatenate_predictions(self, predictions):
        """Concatenate the predictions of each first layer learner.

        This helper is in charge of ensuring the preditions are 2D arrays and
        it will drop one of the probability column when the problem is binary
        since a p(c=1) = 1 - p(c=2).
        """
        X_meta = []
        for est_idx, preds in enumerate(predictions):
            # case where the the estimator returned a 1D array
            if preds.ndim == 1:
                X_meta.append(preds.reshape(-1, 1))
            else:
                # remove one column in case of probabilities with binary case
                if (self.predict_method_[est_idx] == 'predict_proba' and
                        preds.shape[1] == 2):
                    X_meta.append(preds[:, [1]])
                else:
                    X_meta.append(preds)
        return np.concatenate(X_meta, axis=1)

    def fit(self, X, y, sample_weight=None):
        """Fit the estimators.

        Parameters
        ----------
        X : {array-like, sparse matrix}, shape (n_samples, n_features)
            Training vectors, where n_samples is the number of samples and
            n_features is the number of features.

        y : array-like, shape (n_samples,)
            Target values.

        sample_weight : array-like, shape (n_samples,) or None
            Sample weights. If None, then samples are equally weighted.
            Note that this is supported only if all underlying estimators
            support sample weights.

        Returns
        -------
        self : object
        """

        self._validate_meta_estimator()

        if self.estimators is None or len(self.estimators) == 0:
            raise AttributeError(
                "Invalid 'estimators' attribute, 'estimators' should be a list"
                " of (string, estimator) tuples."
            )

        names, estimators_ = zip(*self.estimators)
        self._validate_names(names)

        has_estimator = False
        for est in estimators_:
            if est not in (None, 'drop'):
                has_estimator = True

        if not has_estimator:
            raise ValueError(
                "All estimators are None or 'drop'. At least one is required "
                "to be an estimator."
            )

        if isinstance(self.predict_method, str):
            if self.predict_method != 'auto':
                raise AttributeError(
                    "When 'predict_method' is a string, it should be 'auto'. "
                    " Got {} instead.".format(self.predict_method)
                )
            predict_method = [self.predict_method] * len(estimators_)
        else:
            if len(self.estimators) != len(self.predict_method):
                raise AttributeError(
                    "When 'predict_method' is a list, it should be the same "
                    "length as the list of estimators. Provided {} methods "
                    "for {} estimators."
                    .format(len(self.predict_method), len(self.estimators))
                )
            predict_method = self.predict_method

        # Fit the base estimators on the whole training data. Those
        # base estimators will be used in transform, predict, and
        # predict_proba. They are exposed publicly.
        self.estimators_ = Parallel(n_jobs=self.n_jobs)(
            delayed(_parallel_fit_estimator)(clone(est), X, y, sample_weight)
            for est in estimators_ if est not in (None, 'drop')
        )

        self.named_estimators_ = Bunch()
        est_fitted_idx = 0
        for name_est, org_est in zip(names, estimators_):
            if org_est not in (None, 'drop'):
                self.named_estimators_[name_est] = self.estimators_[
                    est_fitted_idx]
                est_fitted_idx += 1

        # To train the meta-classifier using the most data as possible, we use
        # a cross-validation to obtain the output of the stacked estimators.

        # To ensure that the data provided to each estimator are the same, we
        # need to set the random state of the cv if there is one and we need to
        # take a copy.
        random_state = check_random_state(self.random_state)
        cv = check_cv(self.cv)
        if hasattr(cv, 'random_state'):
            cv.random_state = random_state

        self.predict_method_ = [
            self._method_name(name, est, meth)
            for name, est, meth in zip(names, estimators_, predict_method)
        ]

        predictions = Parallel(n_jobs=self.n_jobs)(
            delayed(cross_val_predict)(clone(est), X, y, cv=deepcopy(cv),
                                       method=meth, n_jobs=self.n_jobs,
                                       verbose=self.verbose)
            for est, meth in zip(estimators_, self.predict_method_)
            if est not in (None, 'drop')
        )

        # Only not None or not 'drop' estimators will be used in transform.
        # Remove the None from the method as well.
        self.predict_method_ = [
            meth for (meth, est) in zip(self.predict_method_, estimators_)
            if est not in (None, 'drop')
        ]

        X_meta = self._concatenate_predictions(predictions)
        self.final_estimator_.fit(X_meta, y)

        return self

    def transform(self, X):
        """Return class labels or probabilities for X for each estimator.

        Parameters
        ----------
        X : {array-like, sparse matrix}, shape (n_samples, n_features)
            Training vectors, where n_samples is the number of samples and
            n_features is the number of features.

        Returns
        -------
        y_preds : ndarray, shape (n_samples, n_estimators)
            Prediction outputs for each estimator.
        """
        check_is_fitted(self, 'estimators_')
        predictions = [
            getattr(est, meth)(X)
            for est, meth in zip(self.estimators_, self.predict_method_)
            if est not in (None, 'drop')
        ]
        return self._concatenate_predictions(predictions)

    @if_delegate_has_method(delegate='final_estimator_')
    def predict(self, X, **predict_params):
        """Predict target for X.

        Parameters
        ----------
        X : {array-like, sparse matrix}, shape (n_samples, n_features)
            Training vectors, where n_samples is the number of samples and
            n_features is the number of features.

        **predict_params : dict of string -> object
            Parameters to the `predict` called by the `final_estimator`. Note
            that this may be used to return uncertainties from some estimators
            with `return_std` or `return_cov`.

        Returns
        -------
        y_pred : ndarray, shape (n_samples,) or (n_samples, n_output)
            Predicted targets.
        """

        check_is_fitted(self, ['estimators_', 'final_estimator_'])
        return self.final_estimator_.predict(
            self.transform(X), **predict_params
        )


class StackingClassifier(_BaseStacking, ClassifierMixin):
    """Stacked of estimators using a final classifier.

    Stacked generalization consists in stacking the output of individual
    estimator and use a classifier to compute the final prediction. Stacking
    allows to use the strength of each individual estimator by using their
    output as input of a final estimator.

    Note that `estimators` are fitted on the full `X` while `final_estimator`
    is trained using cross-validated predictions of the base estimators using
    `cross_val_predict`).

    .. versionadded:: 0.22

    Read more in the :ref:`User Guide <stacking>`.

    Parameters
    ----------
    estimators : list of (string, estimator) tuples
        Base estimators which will be stacked together. An estimator can be set
        to None or 'drop' using `set_params`.

    final_estimator : estimator object, default=None
        A classifier which will be used to combine the base estimators.
        The default classifier is a `LogisticRegression`.

    cv : int, cross-validation generator or an iterable, default=None
        Determines the cross-validation splitting strategy. Possible inputs for
        cv are:

        * None, to use the default 5-fold cross validation,
        * integer, to specify the number of folds in a (Stratified) KFold,
        * An object to be used as a cross-validation generator,
        * An iterable yielding train, test splits.

        For integer/None inputs, if the estimator is a classifier and y is
        either binary or multiclass, `StratifiedKFold` is used. In all other
        cases, `KFold` is used.

        Refer :ref:`User Guide <cross_validation>` for the various
        cross-validation strategies that can be used here.

    predict_method : list of string or 'auto', default='auto'
        Methods called for each base estimator. It can be:

        * if a list of string in which each string is associated to the
          `estimators`,
        * if 'auto', it will try to invoke, for each estimator,
          `predict_proba`, `decision_function` or `predict` in that order.

    n_jobs : int, default=None
        The number of jobs to run in parallel for both `fit` the `estimators`.
        `None` means 1 unless in a `joblib.parallel_backend` context. -1 means
        using all processors. See Glossary for more details.

    random_state : int, RandomState instance or None, optional (default=None)
        If int, random_state is the seed used by the random number generator;
        If RandomState instance, random_state is the random number generator;
        If None, the random number generator is the RandomState instance used
        by `np.random`. Used to set the `cv`.

    Attributes
    ----------
    estimators_ : list of estimator object
        The base estimators fitted.

    named_estimators_ : Bunch object, a dictionary with attribute access
        Attribute to access any fitted sub-estimators by name.

    final_estimator_ : estimator object
        The classifier to stacked the base estimators fitted.

    predict_method_ : list of string
        The method used by each base estimator.

    References
    ----------
    .. [1] Wolpert, David H. "Stacked generalization." Neural networks 5.2
       (1992): 241-259.

    Examples
    --------
    >>> from sklearn.datasets import load_iris
    >>> from sklearn.linear_model import LogisticRegression
    >>> from sklearn.svm import LinearSVC
    >>> from sklearn.ensemble import RandomForestClassifier
    >>> from sklearn.ensemble import StackingClassifier
    >>> X, y = load_iris(return_X_y=True)
    >>> estimators = [
    ...     ('lr', LogisticRegression(solver='lbfgs', multi_class='auto',
    ...                               tol=1e-1)),
    ...     ('svr', LinearSVC(tol=1e-1, random_state=42))
    ... ]
    >>> clf = StackingClassifier(
    ...     estimators=estimators,
    ...     final_estimator=RandomForestClassifier(n_estimators=10,
    ...                                            random_state=42),
    ...     cv=5,
    ... )
    >>> from sklearn.model_selection import train_test_split
    >>> X_train, X_test, y_train, y_test = train_test_split(
    ... X, y, stratify=y, random_state=42
    ... )
    >>> clf.fit(X_train, y_train).score(X_test, y_test)
    0...

    """
    def __init__(self, estimators, final_estimator=None, cv=None,
                 predict_method='auto', n_jobs=None, random_state=None,
                 verbose=0):
        super().__init__(
            estimators=estimators,
            final_estimator=final_estimator,
            cv=cv,
            predict_method=predict_method,
            n_jobs=n_jobs,
            random_state=random_state,
            verbose=verbose
        )

    def _validate_meta_estimator(self):
        super()._validate_final_estimator(
            default=LogisticRegression(random_state=self.random_state)
        )
        if not is_classifier(self.final_estimator_):
            raise AttributeError(
                "'final_estimator' parameter should be a classifier. Got {}"
                .format(self.final_estimator_)
            )

    @if_delegate_has_method(delegate='final_estimator_')
    def predict_proba(self, X):
        """Predict class probabilities for X.

        Parameters
        ----------
        X : {array-like, sparse matrix} of shape (n_samples, n_features)
            Training vectors, where n_samples is the number of samples and
            n_features is the number of features.

        Returns
        -------
        probabilities : ndarray, shape (n_samples, n_classes) or \
            list of arrays n_output
            The class probabilities of the input samples.
        """
        check_is_fitted(self, ['estimators_', 'final_estimator_'])
        return self.final_estimator_.predict_proba(self.transform(X))


class StackingRegressor(_BaseStacking, RegressorMixin):
    """Stacked of estimators using a final regressor.

    Stacked generalization consists in stacking the output of individual
    estimator and use a regressor to compute the final prediction. Stacking
    allows to use the strength of each individual estimator by using their
    output as input of a final estimator.

    Note that `estimators` are fitted on the full `X` while `final_estimator`
    is trained using cross-validated predictions of the base estimators using
    `cross_val_predict`).

    .. versionadded:: 0.22

    Read more in the :ref:`User Guide <stacking>`.

    Parameters
    ----------
    estimators : list of (string, estimator) tuples
        Base estimators which will be stacked together. An estimator can be set
        to None or 'drop' using `set_params`.

    final_estimator : estimator object, default=None
        A regressor which will be used to combine the base estimators.
        The default regressor is a `LinearRegressor`.

    cv : int, cross-validation generator or an iterable, default=None
        Determines the cross-validation splitting strategy. Possible inputs for
        cv are:

        * None, to use the default 5-fold cross validation,
        * integer, to specify the number of folds in a (Stratified) KFold,
        * An object to be used as a cross-validation generator,
        * An iterable yielding train, test splits.

        For integer/None inputs, if the estimator is a classifier and y is
        either binary or multiclass, `StratifiedKFold` is used. In all other
        cases, `KFold` is used.

        Refer :ref:`User Guide <cross_validation>` for the various
        cross-validation strategies that can be used here.

        Refer :ref:`User Guide <cross_validation>` for the various
        cross-validation strategies that can be used here.

    predict_method : list of string or 'auto', default='auto'
        Methods called for each base estimator. It can be:

        * if a list of string in which each string is associated to the
          `estimators`,
        * if 'auto', it will try to invoke, for each estimator,
          `predict_proba`, `decision_function` or `predict` in that order.

    n_jobs : int, default=None
        The number of jobs to run in parallel for both `fit` the `estimators`.
        `None` means 1 unless in a `joblib.parallel_backend` context. -1 means
        using all processors. See Glossary for more details.

    random_state : int, RandomState instance or None, optional (default=None)
        If int, random_state is the seed used by the random number generator;
        If RandomState instance, random_state is the random number generator;
        If None, the random number generator is the RandomState instance used
        by `np.random`. Used to set the `cv`.

    Attributes
    ----------
    estimators_ : list of estimator object
        The base estimators fitted.

    named_estimators_ : Bunch object, a dictionary with attribute access
        Attribute to access any fitted sub-estimators by name.

    final_estimator_ : estimator object
        The regressor to stacked the base estimators fitted.

    predict_method_ : list of string
        The method used by each base estimator.

    References
    ----------
    .. [1] Wolpert, David H. "Stacked generalization." Neural networks 5.2
       (1992): 241-259.

    Examples
    --------
    >>> from sklearn.datasets import load_diabetes
    >>> from sklearn.linear_model import LinearRegression
    >>> from sklearn.svm import LinearSVR
    >>> from sklearn.ensemble import RandomForestRegressor
    >>> from sklearn.ensemble import StackingRegressor
    >>> estimators = [
    ...     ('lr', LinearRegression()),
    ...     ('svr', LinearSVR(tol=1e-1, random_state=42))
    ... ]
    >>> reg = StackingRegressor(
    ...     estimators=estimators,
    ...     final_estimator=RandomForestRegressor(n_estimators=10,
    ...                                           random_state=42),
    ...     cv=5
    ... )
    >>> from sklearn.model_selection import train_test_split
    >>> X_train, X_test, y_train, y_test = train_test_split(
    ... *load_diabetes(return_X_y=True), random_state=42
    ... )
    >>> reg.fit(X_train, y_train).score(X_test, y_test)
    0...

    """
    def __init__(self, estimators, final_estimator=None, cv=None,
                 predict_method='auto', n_jobs=None, random_state=None,
                 verbose=0):
        super().__init__(
            estimators=estimators,
            final_estimator=final_estimator,
            cv=cv,
            predict_method=predict_method,
            n_jobs=n_jobs,
            random_state=random_state,
            verbose=verbose
        )

    def _validate_meta_estimator(self):
        super()._validate_final_estimator(
            default=LinearRegression()
        )
        if not is_regressor(self.final_estimator_):
            raise AttributeError(
                "'final_estimator' parameter should be a regressor. Got {}"
                .format(self.final_estimator_)
            )
