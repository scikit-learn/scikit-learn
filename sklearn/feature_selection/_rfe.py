# Authors: Alexandre Gramfort <alexandre.gramfort@inria.fr>
#          Vincent Michel <vincent.michel@inria.fr>
#          Gilles Louppe <g.louppe@gmail.com>
#
# License: BSD 3 clause

"""Recursive feature elimination for feature ranking"""

import numpy as np
import numbers
from joblib import Parallel, delayed, effective_n_jobs
import matplotlib.pyplot as plt

from ..utils.metaestimators import if_delegate_has_method
from ..utils.metaestimators import _safe_split
from ..utils.validation import check_is_fitted
from ..utils.validation import _deprecate_positional_args
from ..base import BaseEstimator
from ..base import MetaEstimatorMixin
from ..base import clone
from ..base import is_classifier
from ..model_selection import check_cv
from ..model_selection._validation import _score
from ..metrics import check_scoring
from ._base import SelectorMixin
from ._base import _get_feature_importances


def _rfe_single_fit(rfe, estimator, X, y, train, test, scorer):
    """
    Return the score for a fit across one fold.
    """
    X_train, y_train = _safe_split(estimator, X, y, train)
    X_test, y_test = _safe_split(estimator, X, y, test, train)
    return rfe._fit(
        X_train, y_train, lambda estimator, features:
        _score(estimator, X_test[:, features], y_test, scorer)).scores_


def std_interval(array):
    '''
    Return the standard diviation intercval of an arrary.
    '''
    mean = np.mean(array)
    std = np.std(array)
    return(mean-std, mean+std)


class RFE(SelectorMixin, MetaEstimatorMixin, BaseEstimator):
    """Feature ranking with recursive feature elimination.

    Given an external estimator that assigns weights to features (e.g., the
    coefficients of a linear model), the goal of recursive feature elimination
    (RFE) is to select features by recursively considering smaller and smaller
    sets of features. First, the estimator is trained on the initial set of
    features and the importance of each feature is obtained either through
    any specific attribute or callable.
    Then, the least important features are pruned from current set of features.
    That procedure is recursively repeated on the pruned set until the desired
    number of features to select is eventually reached.

    Read more in the :ref:`User Guide <rfe>`.

    Parameters
    ----------
    estimator : ``Estimator`` instance
        A supervised learning estimator with a ``fit`` method that provides
        information about feature importance
        (e.g. `coef_`, `feature_importances_`).

    n_features_to_select : int or float, default=None
        The number of features to select. If `None`, half of the features are
        selected. If integer, the parameter is the absolute number of features
        to select. If float between 0 and 1, it is the fraction of features to
        select.

    step : int or float, default=1
        If greater than or equal to 1, then ``step`` corresponds to the
        (integer) number of features to remove at each iteration.
        If within (0.0, 1.0), then ``step`` corresponds to the percentage
        (rounded down) of features to remove at each iteration.

    verbose : int, default=0
        Controls verbosity of output.

    importance_getter : str or callable, default='auto'
        If 'auto', uses the feature importance either through a `coef_`
        or `feature_importances_` attributes of estimator.

        Also accepts a string that specifies an attribute name/path
        for extracting feature importance (implemented with `attrgetter`).
        For example, give `regressor_.coef_` in case of
        :class:`~sklearn.compose.TransformedTargetRegressor`  or
        `named_steps.clf.feature_importances_` in case of
        class:`~sklearn.pipeline.Pipeline` with its last step named `clf`.

        If `callable`, overrides the default feature importance getter.
        The callable is passed with the fitted estimator and it should
        return importance for each feature.

        .. versionadded:: 0.24

    Attributes
    ----------
    estimator_ : ``Estimator`` instance
        The fitted estimator used to select features.

    n_features_ : int
        The number of selected features.

    ranking_ : ndarray of shape (n_features,)
        The feature ranking, such that ``ranking_[i]`` corresponds to the
        ranking position of the i-th feature. Selected (i.e., estimated
        best) features are assigned rank 1.

    support_ : ndarray of shape (n_features,)
        The mask of selected features.

    Examples
    --------
    The following example shows how to retrieve the 5 most informative
    features in the Friedman #1 dataset.

    >>> from sklearn.datasets import make_friedman1
    >>> from sklearn.feature_selection import RFE
    >>> from sklearn.svm import SVR
    >>> X, y = make_friedman1(n_samples=50, n_features=10, random_state=0)
    >>> estimator = SVR(kernel="linear")
    >>> selector = RFE(estimator, n_features_to_select=5, step=1)
    >>> selector = selector.fit(X, y)
    >>> selector.support_
    array([ True,  True,  True,  True,  True, False, False, False, False,
           False])
    >>> selector.ranking_
    array([1, 1, 1, 1, 1, 6, 4, 3, 2, 5])

    Notes
    -----
    Allows NaN/Inf in the input if the underlying estimator does as well.

    See also
    --------
    RFECV : Recursive feature elimination with built-in cross-validated
        selection of the best number of features

    References
    ----------

    .. [1] Guyon, I., Weston, J., Barnhill, S., & Vapnik, V., "Gene selection
           for cancer classification using support vector machines",
           Mach. Learn., 46(1-3), 389--422, 2002.
    """
    @_deprecate_positional_args
    def __init__(self, estimator, *, n_features_to_select=None, step=1,
                 verbose=0, importance_getter='auto'):
        self.estimator = estimator
        self.n_features_to_select = n_features_to_select
        self.step = step
        self.importance_getter = importance_getter
        self.verbose = verbose

    @property
    def _estimator_type(self):
        return self.estimator._estimator_type

    @property
    def classes_(self):
        return self.estimator_.classes_

    def fit(self, X, y):
        """Fit the RFE model and then the underlying estimator on the selected
           features.

        Parameters
        ----------
        X : {array-like, sparse matrix} of shape (n_samples, n_features)
            The training input samples.

        y : array-like of shape (n_samples,)
            The target values.
        """
        return self._fit(X, y)

    def _fit(self, X, y, step_score=None):
        # Parameter step_score controls the calculation of self.scores_
        # step_score is not exposed to users
        # and is used when implementing RFECV
        # self.scores_ will not be calculated when calling _fit through fit

        tags = self._get_tags()
        X, y = self._validate_data(
            X, y, accept_sparse="csc",
            ensure_min_features=2,
            force_all_finite=not tags.get('allow_nan', True),
            multi_output=True
        )
        error_msg = ("n_features_to_select must be either None, a "
                     "positive integer representing the absolute "
                     "number of features or a float in (0.0, 1.0] "
                     "representing a percentage of features to "
                     f"select. Got {self.n_features_to_select}")

        # Initialization
        n_features = X.shape[1]
        if self.n_features_to_select is None:
            n_features_to_select = n_features // 2
        elif self.n_features_to_select < 0:
            raise ValueError(error_msg)
        elif isinstance(self.n_features_to_select, numbers.Integral):  # int
            n_features_to_select = self.n_features_to_select
        elif self.n_features_to_select > 1.0:  # float > 1
            raise ValueError(error_msg)
        else:  # float
            n_features_to_select = int(n_features * self.n_features_to_select)

        if 0.0 < self.step < 1.0:
            step = int(max(1, self.step * n_features))
        else:
            step = int(self.step)
        if step <= 0:
            raise ValueError("Step must be >0")

        support_ = np.ones(n_features, dtype=bool)
        ranking_ = np.ones(n_features, dtype=int)

        if step_score:
            self.scores_ = []

        # Elimination
        while np.sum(support_) > n_features_to_select:
            # Remaining features
            features = np.arange(n_features)[support_]

            # Rank the remaining features
            estimator = clone(self.estimator)
            if self.verbose > 0:
                print("Fitting estimator with %d features." % np.sum(support_))

            estimator.fit(X[:, features], y)

            # Get importance and rank them
            importances = _get_feature_importances(
                estimator, self.importance_getter, transform_func="square",
            )
            ranks = np.argsort(importances)

            # for sparse case ranks is matrix
            ranks = np.ravel(ranks)

            # Eliminate the worse features
            threshold = min(step, np.sum(support_) - n_features_to_select)

            # Compute step score on the previous selection iteration
            # because 'estimator' must use features
            # that have not been eliminated yet
            if step_score:
                self.scores_.append(step_score(estimator, features))
            support_[features[ranks][:threshold]] = False
            ranking_[np.logical_not(support_)] += 1

        # Set final attributes
        features = np.arange(n_features)[support_]
        self.estimator_ = clone(self.estimator)
        self.estimator_.fit(X[:, features], y)

        # Compute step score when only n_features_to_select features left
        if step_score:
            self.scores_.append(step_score(self.estimator_, features))
        self.n_features_ = support_.sum()
        self.support_ = support_
        self.ranking_ = ranking_

        return self

    @if_delegate_has_method(delegate='estimator')
    def predict(self, X):
        """Reduce X to the selected features and then predict using the
           underlying estimator.

        Parameters
        ----------
        X : array of shape [n_samples, n_features]
            The input samples.

        Returns
        -------
        y : array of shape [n_samples]
            The predicted target values.
        """
        check_is_fitted(self)
        return self.estimator_.predict(self.transform(X))

    @if_delegate_has_method(delegate='estimator')
    def score(self, X, y):
        """Reduce X to the selected features and then return the score of the
           underlying estimator.

        Parameters
        ----------
        X : array of shape [n_samples, n_features]
            The input samples.

        y : array of shape [n_samples]
            The target values.
        """
        check_is_fitted(self)
        return self.estimator_.score(self.transform(X), y)

    def _get_support_mask(self):
        check_is_fitted(self)
        return self.support_

    @if_delegate_has_method(delegate='estimator')
    def decision_function(self, X):
        """Compute the decision function of ``X``.

        Parameters
        ----------
        X : {array-like or sparse matrix} of shape (n_samples, n_features)
            The input samples. Internally, it will be converted to
            ``dtype=np.float32`` and if a sparse matrix is provided
            to a sparse ``csr_matrix``.

        Returns
        -------
        score : array, shape = [n_samples, n_classes] or [n_samples]
            The decision function of the input samples. The order of the
            classes corresponds to that in the attribute :term:`classes_`.
            Regression and binary classification produce an array of shape
            [n_samples].
        """
        check_is_fitted(self)
        return self.estimator_.decision_function(self.transform(X))

    @if_delegate_has_method(delegate='estimator')
    def predict_proba(self, X):
        """Predict class probabilities for X.

        Parameters
        ----------
        X : {array-like or sparse matrix} of shape (n_samples, n_features)
            The input samples. Internally, it will be converted to
            ``dtype=np.float32`` and if a sparse matrix is provided
            to a sparse ``csr_matrix``.

        Returns
        -------
        p : array of shape (n_samples, n_classes)
            The class probabilities of the input samples. The order of the
            classes corresponds to that in the attribute :term:`classes_`.
        """
        check_is_fitted(self)
        return self.estimator_.predict_proba(self.transform(X))

    @if_delegate_has_method(delegate='estimator')
    def predict_log_proba(self, X):
        """Predict class log-probabilities for X.

        Parameters
        ----------
        X : array of shape [n_samples, n_features]
            The input samples.

        Returns
        -------
        p : array of shape (n_samples, n_classes)
            The class log-probabilities of the input samples. The order of the
            classes corresponds to that in the attribute :term:`classes_`.
        """
        check_is_fitted(self)
        return self.estimator_.predict_log_proba(self.transform(X))

    def _more_tags(self):
        estimator_tags = self.estimator._get_tags()
        return {'poor_score': True,
                'allow_nan': estimator_tags.get('allow_nan', True),
                'requires_y': True,
                }


class RFECV(RFE):
    """Feature ranking with recursive feature elimination and cross-validated
    selection of the best number of features.

    See glossary entry for :term:`cross-validation estimator`.

    Read more in the :ref:`User Guide <rfe>`.

    Parameters
    ----------
    estimator : ``Estimator`` instance
        A supervised learning estimator with a ``fit`` method that provides
        information about feature importance either through a ``coef_``
        attribute or through a ``feature_importances_`` attribute.

    step : int or float, default=1
        If greater than or equal to 1, then ``step`` corresponds to the
        (integer) number of features to remove at each iteration.
        If within (0.0, 1.0), then ``step`` corresponds to the percentage
        (rounded down) of features to remove at each iteration.
        Note that the last iteration may remove fewer than ``step`` features in
        order to reach ``min_features_to_select``.

    min_features_to_select : int, default=1
        The minimum number of features to be selected. This number of features
        will always be scored, even if the difference between the original
        feature count and ``min_features_to_select`` isn't divisible by
        ``step``.

        .. versionadded:: 0.20

    cv : int, cross-validation generator or an iterable, default=None
        Determines the cross-validation splitting strategy.
        Possible inputs for cv are:

        - None, to use the default 5-fold cross-validation,
        - integer, to specify the number of folds.
        - :term:`CV splitter`,
        - An iterable yielding (train, test) splits as arrays of indices.

        For integer/None inputs, if ``y`` is binary or multiclass,
        :class:`~sklearn.model_selection.StratifiedKFold` is used. If the
        estimator is a classifier or if ``y`` is neither binary nor multiclass,
        :class:`~sklearn.model_selection.KFold` is used.

        Refer :ref:`User Guide <cross_validation>` for the various
        cross-validation strategies that can be used here.

        .. versionchanged:: 0.22
            ``cv`` default value of None changed from 3-fold to 5-fold.

    scoring : string, callable or None, default=None
        A string (see model evaluation documentation) or
        a scorer callable object / function with signature
        ``scorer(estimator, X, y)``.

    verbose : int, default=0
        Controls verbosity of output.

    n_jobs : int or None, default=None
        Number of cores to run in parallel while fitting across folds.
        ``None`` means 1 unless in a :obj:`joblib.parallel_backend` context.
        ``-1`` means using all processors. See :term:`Glossary <n_jobs>`
        for more details.

        .. versionadded:: 0.18

    importance_getter : str or callable, default='auto'
        If 'auto', uses the feature importance either through a `coef_`
        or `feature_importances_` attributes of estimator.

        Also accepts a string that specifies an attribute name/path
        for extracting feature importance.
        For example, give `regressor_.coef_` in case of
        :class:`~sklearn.compose.TransformedTargetRegressor`  or
        `named_steps.clf.feature_importances_` in case of
        :class:`~sklearn.pipeline.Pipeline` with its last step named `clf`.

        If `callable`, overrides the default feature importance getter.
        The callable is passed with the fitted estimator and it should
        return importance for each feature.

        .. versionadded:: 0.24

    active_feature_size_selection: bool, default='False'
        if True, rankings will generate for every feature and will enable
        selected_feature_size parameter in the adjust_feature_size method.
        if False, normal rfecv will be performed.

    Attributes
    ----------
    estimator_ : ``Estimator`` instance
        The fitted estimator used to select features.

    grid_scores_ : ndarray of shape (n_subsets_of_features,)
        The cross-validation scores such that
        ``grid_scores_[i]`` corresponds to
        the CV score of the i-th subset of features.

    n_features_ : int
        The number of selected features with cross-validation.

    ranking_ : narray of shape (n_features,)
        The feature ranking, such that `ranking_[i]`
        corresponds to the ranking
        position of the i-th feature.
        Selected (i.e., estimated best)
        features are assigned rank 1.

    support_ : ndarray of shape (n_features,)
        The mask of selected features.

    Notes
    -----
    The size of ``grid_scores_`` is equal to
    ``ceil((n_features - min_features_to_select) / step) + 1``,
    where step is the number of features removed at each iteration.

    Allows NaN/Inf in the input if the underlying estimator does as well.

    Examples
    --------
    The following example shows how to retrieve the a-priori not known 5
    informative features in the Friedman #1 dataset.

    >>> from sklearn.datasets import make_friedman1
    >>> from sklearn.feature_selection import RFECV
    >>> from sklearn.svm import SVR
    >>> X, y = make_friedman1(n_samples=50, n_features=10, random_state=0)
    >>> estimator = SVR(kernel="linear")
    >>> selector = RFECV(estimator, step=1, cv=5)
    >>> selector = selector.fit(X, y)
    >>> selector.support_
    array([ True,  True,  True,  True,  True, False, False, False, False,
           False])
    >>> selector.ranking_
    array([1, 1, 1, 1, 1, 6, 4, 3, 2, 5])

    See also
    --------
    RFE : Recursive feature elimination

    References
    ----------

    .. [1] Guyon, I., Weston, J., Barnhill, S., & Vapnik, V., "Gene selection
           for cancer classification using support vector machines",
           Mach. Learn., 46(1-3), 389--422, 2002.
    """
    @_deprecate_positional_args
    def __init__(self, estimator, *, step=1, min_features_to_select=1,
                 cv=None, scoring=None, verbose=0, n_jobs=None,
                 importance_getter='auto',
                 active_feature_size_selection=False):
        self.estimator = estimator
        self.step = step
        self.importance_getter = importance_getter
        self.cv = cv
        self.scoring = scoring
        self.verbose = verbose
        self.n_jobs = n_jobs
        self.min_features_to_select = min_features_to_select
        self.active_feature_size_selection = active_feature_size_selection

    def fit(self, X, y, groups=None):
        """Fit the RFE model and automatically tune the number of selected
           features.

        Parameters
        ----------
        X : {array-like, sparse matrix} of shape (n_samples, n_features)
            Training vector, where `n_samples` is the number of samples and
            `n_features` is the total number of features.

        y : array-like of shape (n_samples,)
            Target values (integers for classification, real numbers for
            regression).

        groups : array-like of shape (n_samples,) or None, default=None
            Group labels for the samples used while splitting the dataset into
            train/test set. Only used in conjunction with a "Group" :term:`cv`
            instance (e.g., :class:`~sklearn.model_selection.GroupKFold`).

            .. versionadded:: 0.20
        """
        tags = self._get_tags()
        X, y = self._validate_data(
            X, y, accept_sparse="csr", ensure_min_features=2,
            force_all_finite=not tags.get('allow_nan', True),
            multi_output=True
        )

        # Initialization
        cv = check_cv(self.cv, y, classifier=is_classifier(self.estimator))
        scorer = check_scoring(self.estimator, scoring=self.scoring)
        n_features = X.shape[1]

        if 0.0 < self.step < 1.0:
            step = int(max(1, self.step * n_features))
        else:
            step = int(self.step)
        if step <= 0:
            raise ValueError("Step must be >0")

        # Build an RFE object, which will evaluate and score each possible
        # feature count, down to self.min_features_to_select
        rfe = RFE(estimator=self.estimator,
                  n_features_to_select=self.min_features_to_select,
                  importance_getter=self.importance_getter,
                  step=self.step, verbose=self.verbose)

        # Determine the number of subsets of features by fitting across
        # the train folds and choosing the "features_to_select" parameter
        # that gives the least averaged error across all folds.

        # Note that joblib raises a non-picklable error for bound methods
        # even if n_jobs is set to 1 with the default multiprocessing
        # backend.
        # This branching is done so that to
        # make sure that user code that sets n_jobs to 1
        # and provides bound methods as scorers is not broken with the
        # addition of n_jobs parameter in version 0.18.

        if effective_n_jobs(self.n_jobs) == 1:
            parallel, func = list, _rfe_single_fit
        else:
            parallel = Parallel(n_jobs=self.n_jobs)
            func = delayed(_rfe_single_fit)

        scores = parallel(
            func(rfe, self.estimator, X, y, train, test, scorer)
            for train, test in cv.split(X, y, groups))

        scores = np.sum(scores, axis=0)
        scores_rev = scores[::-1]
        argmax_idx = len(scores) - np.argmax(scores_rev) - 1

        # Set the minimum number of features to 1
        # in order to record the ranking for every feature
        if self.active_feature_size_selection:
            n_features_to_select = 1

        # If active_feature_size_selection is False
        # Nothing will be changed
        else:
            n_features_to_select = max(
                n_features - (argmax_idx * step),
                self.min_features_to_select)

        # Re-execute an elimination with best_k over the whole set
        rfe = RFE(estimator=self.estimator,
                  n_features_to_select=n_features_to_select, step=self.step,
                  importance_getter=self.importance_getter,
                  verbose=self.verbose)

        rfe.fit(X, y)

        # Set final attributes
        self.support_ = rfe.support_
        self.n_features_ = rfe.n_features_
        self.ranking_ = rfe.ranking_
        self.estimator_ = clone(self.estimator)
        self.estimator_.fit(self.transform(X), y)

        # Fixing a normalization error, n is equal to get_n_splits(X, y) - 1
        # here, the scores are normalized by get_n_splits(X, y)
        self.grid_scores_ = scores[::-1] / cv.get_n_splits(X, y, groups)

        # Set attributes to optimal feature size
        # if active_feature_size_selection is true
        if self.active_feature_size_selection:
            self.ranking_original = self.ranking_
            optimal_feature_size = list(self.grid_scores_).index(
                np.max(self.grid_scores_))+1
            self.ranking_ = np.array(
                [i-optimal_feature_size+1 if i > optimal_feature_size
                 else 1 for i in self.ranking_original])
            self.support_ = np.array(
                [True if i == 1 else False for i in self.ranking_])
            self.n_features_ = list(self.ranking_).count(1)

        return self

    def adjust_feature_size(self, selected_feature_size=None):
        '''Adjust the feature size customizablly with visualization.
        (from 1 to the full feature size)
        (With updated ranking_, support_, n_features attributes)

        Parameters
        ----------        
        selected_feature_size : int, (default=None)
            The chosen number of features selected. If the number is greater 
            than the full feature size or if the variable is set to None, only 
            the optimal feature size will be shown.

        Examples
        --------
        The following example shows how to chose the custom feature size with 
        the adjust_feature_size function (Minimal feature size first)

        # Import library
        from sklearn.ensemble import RandomForestClassifier
        from sklearn.datasets import make_classification

        # Make data
        X, y = make_classification(n_samples=50,
                                   n_features=30,
                                   n_informative=5,
                                   n_redundant=0,
                                   n_repeated=0,
                                   n_classes=2,
                                   random_state=0,
                                   shuffle=False)

        # Build a random forest model and 
        # perform the RFE with 5-fold cross validation
        clf = RandomForestClassifier(random_state=0)
        rfecv = RFECV(estimator=clf, step=1, cv=StratifiedKFold(5),
                      scoring='roc_auc',n_jobs=-1,
                      active_feature_size_selection=True)
        rfecv.fit(X, y)

        # Adjust the feature size to maximize your need
        rfecv.adjust_feature_size(selected_feature_size=26)
        rfecv.adjust_feature_size(selected_feature_size=15)

        # Access the updated attitbute
        print(rfecv.ranking_)'''

        # calculate mean score, standard deviation,
        # score intervals and optimal feature_size
        optimal_feature_size = list(self.grid_scores_).index(
            np.max(self.grid_scores_))+1
        score_lower, score_upper = zip(
            *[std_interval(row) for row in list(zip(*self.scores))[::-1]])
        n_features_selected = range(1, len(self.grid_scores_) + 1)

        # Visualize
        plt.figure(figsize=(9, 6))
        plt.plot(n_features_selected, self.grid_scores_, color='#14213d',
                 label=r'Mean score',
                 lw=2, alpha=.8)
        plt.fill_between(n_features_selected, score_lower,
                         score_upper, color='#90a8c3', alpha=.2,
                         label=r'$\pm$ 1 std. dev.')

        if selected_feature_size and self.active_feature_size_selection:
            # Requires self.active_feature_size_selection to be True
            # change attribute
            self.ranking_ = np.array(
                [i-selected_feature_size+1 if i > selected_feature_size
                 else 1 for i in self.ranking_original])
            self.support_ = np.array(
                [True if i == 1 else False for i in self.ranking_])
            self.n_features_ = list(self.ranking_).count(1)
            # plot
            selected_size_score = self.grid_scores_[selected_feature_size-1]
            plt.axvline(x=selected_feature_size, linestyle='--',
                        lw=2, color='#b2182b',
                        label='Selected feature size \n\
                        (score = %0.3f) \n(feature size = %i)'
                        % (selected_size_score,
                            selected_feature_size), alpha=.8)
            plt.axvline(x=optimal_feature_size, linestyle='--',
                        lw=2, color='#d7b9d5',
                        label='Optimal feature size \n\
                        (score = %0.3f) \n(feature size = %i)'
                        % (np.max(self.grid_scores_),
                            optimal_feature_size), alpha=.8)

        else:
            if not self.active_feature_size_selection \
                    and selected_feature_size:
                print('In order to select feature size, \
                    please set active_feature_size_selection to True')
            plt.axvline(x=optimal_feature_size, linestyle='--',
                        lw=2, color='#b2182b',
                        label='Optimal feature size \n\
                        (score = %0.3f) \n(feature size = %i)'
                        % (np.max(self.grid_scores_),
                            optimal_feature_size), alpha=.8)
        plt.xlabel("Number of features selected")
        plt.ylabel("Cross validation score (nb of correct classifications)")
        plt.title(f'Recursive Feature Elimination with Cross Validation')
        plt.legend(loc="lower right", labelspacing=1.3)
        plt.show()
