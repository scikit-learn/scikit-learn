# Authors: Alexandre Gramfort <alexandre.gramfort@inria.fr>
#          Vincent Michel <vincent.michel@inria.fr>
#          Gilles Louppe <g.louppe@gmail.com>
#          Leandro Hermida <hermidal@cs.umd.edu>
#
# License: BSD 3 clause

"""Recursive feature elimination for feature ranking"""

import numpy as np
from ..utils import check_X_y, safe_sqr
from ..utils.metaestimators import if_delegate_has_method
from ..utils.metaestimators import _safe_split
from ..utils.validation import check_is_fitted
from ..base import BaseEstimator
from ..base import MetaEstimatorMixin
from ..base import clone
from ..base import is_classifier
from ..utils._joblib import Parallel, delayed, effective_n_jobs
from ..model_selection import check_cv
from ..model_selection._validation import _score
from ..metrics.scorer import check_scoring
from .base import SelectorMixin


def _rfe_single_fit(rfe, estimator, X, y, train, test, scorer):
    """
    Return the score for a fit across one fold.
    """
    X_train, y_train = _safe_split(estimator, X, y, train)
    X_test, y_test = _safe_split(estimator, X, y, test, train)
    rfe._fit(
        X_train, y_train, lambda estimator, features:
        _score(estimator, X_test[:, features], y_test, scorer))
    return rfe.scores_, rfe.n_remaining_feature_steps_


class RFE(BaseEstimator, MetaEstimatorMixin, SelectorMixin):
    """Feature ranking with recursive feature elimination.

    Given an external estimator that assigns weights to features (e.g., the
    coefficients of a linear model), the goal of recursive feature elimination
    (RFE) is to select features by recursively considering smaller and smaller
    sets of features. First, the estimator is trained on the initial set of
    features and the importance of each feature is obtained either through a
    ``coef_`` attribute or through a ``feature_importances_`` attribute.
    Then, the least important features are pruned from current set of features.
    That procedure is recursively repeated on the pruned set until the desired
    number of features to select is eventually reached.

    Read more in the :ref:`User Guide <rfe>`.

    Parameters
    ----------
    estimator : object
        A supervised learning estimator with a ``fit`` method that provides
        information about feature importance either through a ``coef_``
        attribute or through a ``feature_importances_`` attribute.

    n_features_to_select : int or None (default=None)
        The number of features to select. If `None`, half of the features
        are selected.

    step : int or float, optional (default=1)
        If greater than or equal to 1, then ``step`` corresponds to the
        (integer) number of features to remove at each iteration.
        If within (0.0, 1.0), then ``step`` corresponds to the percentage
        (rounded down) of features to remove at each iteration.

    step_change_at : int or float or None, optional (default=None)
        If specified ``step`` int or float equates to > 1 feature, then the
        threshold number of remaining features when to change to
        ``changed_step``. If ``step_change_at`` within (0.0, 1.0), then
        threshold number of remaining features calculated as fraction of
        original number of features.

    changed_step : int or float, optional (default=1)
        Step size to change to starting at ``step_change_at`` if specified. If
        greater than or equal to 1, then ``changed_step`` corresponds to the
        (integer) number of features to remove at each iteration. If within
        (0.0, 1.0), then ``changed_step`` corresponds to the percentage
        (rounded down) of features to remove at each iteration.

    reducing_step : boolean, optional (default=False)
        If true and ``step`` or ``changed_step`` is a float, the number of
        features removed is calculated as a fraction of the remaining features
        in that iteration. If false, the number of features removed is constant
        (a fraction of the original number of features) across iterations.

    verbose : int, (default=0)
        Controls verbosity of output.

    Attributes
    ----------
    n_features_ : int
        The number of selected features.

    support_ : array of shape [n_features]
        The mask of selected features.

    ranking_ : array of shape [n_features]
        The feature ranking, such that ``ranking_[i]`` corresponds to the
        ranking position of the i-th feature. Selected (i.e., estimated
        best) features are assigned rank 1.

    estimator_ : object
        The external estimator fit on the reduced dataset.

    Examples
    --------
    The following example shows how to retrieve the 5 right informative
    features in the Friedman #1 dataset.

    >>> from sklearn.datasets import make_friedman1
    >>> from sklearn.feature_selection import RFE
    >>> from sklearn.svm import SVR
    >>> X, y = make_friedman1(n_samples=50, n_features=10, random_state=0)
    >>> estimator = SVR(kernel="linear")
    >>> selector = RFE(estimator, 5, step=1)
    >>> selector = selector.fit(X, y)
    >>> selector.support_ # doctest: +NORMALIZE_WHITESPACE
    array([ True,  True,  True,  True,  True, False, False, False, False,
           False])
    >>> selector.ranking_
    array([1, 1, 1, 1, 1, 6, 4, 3, 2, 5])

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
    def __init__(self, estimator, n_features_to_select=None, step=1,
                 step_change_at=None, changed_step=1, reducing_step=False,
                 verbose=0):
        self.estimator = estimator
        self.n_features_to_select = n_features_to_select
        self.step = step
        self.step_change_at = step_change_at
        self.changed_step = changed_step
        self.reducing_step = reducing_step
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
        X : {array-like, sparse matrix}, shape = [n_samples, n_features]
            The training input samples.

        y : array-like, shape = [n_samples]
            The target values.
        """
        return self._fit(X, y)

    def _fit(self, X, y, step_score=None):
        # Parameter step_score controls the calculation of self.scores_
        # step_score is not exposed to users
        # and is used when implementing RFECV
        # self.scores_ will not be calculated when calling _fit through fit

        X, y = check_X_y(X, y, "csc", ensure_min_features=2)
        # Initialization
        n_features = X.shape[1]
        if self.n_features_to_select is None:
            n_features_to_select = n_features // 2
        else:
            n_features_to_select = self.n_features_to_select

        if self.step >= 1.0:
            step = int(self.step)
        elif 0.0 < self.step < 1.0 and not self.reducing_step:
            step = int(max(1, self.step * n_features))
        elif self.step <= 0:
            raise ValueError("step must be > 0")

        if self.step_change_at is not None:
            if self.step_change_at >= 1.0:
                step_change_at = int(self.step_change_at)
            elif 0.0 < self.step_change_at < 1.0:
                step_change_at = int(max(1, self.step_change_at * n_features))
            if not n_features_to_select < step_change_at < n_features:
                raise ValueError("step_change_at must be greater than "
                                 "n_features_to_select and less than initial "
                                 "number of features")
            if self.changed_step >= 1.0:
                changed_step = int(self.changed_step)
            elif 0.0 < self.changed_step < 1.0 and not self.reducing_step:
                changed_step = int(max(1, self.changed_step * step_change_at))
            elif self.changed_step <= 0:
                raise ValueError("changed_step must be > 0")

        support_ = np.ones(n_features, dtype=np.bool)
        ranking_ = np.ones(n_features, dtype=np.int)

        if step_score:
            self.scores_ = []
            self.n_remaining_feature_steps_ = []

        # Elimination
        n_remaining_features = n_features
        while n_remaining_features > n_features_to_select:
            # Remaining features
            features = np.arange(n_features)[support_]

            # Rank the remaining features
            estimator = clone(self.estimator)
            if self.verbose > 0:
                print("Fitting estimator with %d features."
                      % n_remaining_features)

            estimator.fit(X[:, features], y)

            # Get coefs
            if hasattr(estimator, 'coef_'):
                coefs = estimator.coef_
            else:
                coefs = getattr(estimator, 'feature_importances_', None)
            if coefs is None:
                raise RuntimeError('The classifier does not expose '
                                   '"coef_" or "feature_importances_" '
                                   'attributes')

            # Get ranks
            if coefs.ndim > 1:
                ranks = np.argsort(safe_sqr(coefs).sum(axis=0))
            else:
                ranks = np.argsort(safe_sqr(coefs))

            # for sparse case ranks is matrix
            ranks = np.ravel(ranks)

            # Adjust step using special parameters if specified
            if self.step_change_at is not None:
                if n_remaining_features > step_change_at:
                    if 0.0 < self.step < 1.0 and self.reducing_step:
                        step = int(max(1, min(
                            n_remaining_features - step_change_at,
                            self.step * n_remaining_features
                        )))
                    else:
                        step = min(
                            n_remaining_features - step_change_at,
                            step
                        )
                elif 0.0 < self.changed_step < 1.0 and self.reducing_step:
                    step = int(max(1, min(
                        n_remaining_features - n_features_to_select,
                        self.changed_step * n_remaining_features
                    )))
                else:
                    step = changed_step
            elif 0.0 < self.step < 1.0 and self.reducing_step:
                step = int(max(1, min(
                    n_remaining_features - n_features_to_select,
                    self.step * n_remaining_features
                )))

            # Eliminate the worse features
            threshold = min(step, n_remaining_features - n_features_to_select)

            # Compute step score on the previous selection iteration
            # because 'estimator' must use features
            # that have not been eliminated yet
            if step_score:
                self.scores_.append(step_score(estimator, features))
                self.n_remaining_feature_steps_.append(n_remaining_features)
            support_[features[ranks][:threshold]] = False
            ranking_[np.logical_not(support_)] += 1
            n_remaining_features -= threshold

        # Set final attributes
        features = np.arange(n_features)[support_]
        self.estimator_ = clone(self.estimator)
        self.estimator_.fit(X[:, features], y)

        # Compute step score when only n_features_to_select features left
        if step_score:
            self.scores_.append(step_score(self.estimator_, features))
            self.n_remaining_feature_steps_.append(n_remaining_features)
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
        check_is_fitted(self, 'estimator_')
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
        check_is_fitted(self, 'estimator_')
        return self.estimator_.score(self.transform(X), y)

    def _get_support_mask(self):
        check_is_fitted(self, 'support_')
        return self.support_

    @if_delegate_has_method(delegate='estimator')
    def decision_function(self, X):
        """Compute the decision function of ``X``.

        Parameters
        ----------
        X : array-like or sparse matrix, shape = [n_samples, n_features]
            The input samples. Internally, it will be converted to
            ``dtype=np.float32`` and if a sparse matrix is provided
            to a sparse ``csr_matrix``.

        Returns
        -------
        score : array, shape = [n_samples, n_classes] or [n_samples]
            The decision function of the input samples. The order of the
            classes corresponds to that in the attribute `classes_`.
            Regression and binary classification produce an array of shape
            [n_samples].
        """
        check_is_fitted(self, 'estimator_')
        return self.estimator_.decision_function(self.transform(X))

    @if_delegate_has_method(delegate='estimator')
    def predict_proba(self, X):
        """Predict class probabilities for X.

        Parameters
        ----------
        X : array-like or sparse matrix, shape = [n_samples, n_features]
            The input samples. Internally, it will be converted to
            ``dtype=np.float32`` and if a sparse matrix is provided
            to a sparse ``csr_matrix``.

        Returns
        -------
        p : array of shape = [n_samples, n_classes]
            The class probabilities of the input samples. The order of the
            classes corresponds to that in the attribute `classes_`.
        """
        check_is_fitted(self, 'estimator_')
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
        p : array of shape = [n_samples, n_classes]
            The class log-probabilities of the input samples. The order of the
            classes corresponds to that in the attribute `classes_`.
        """
        check_is_fitted(self, 'estimator_')
        return self.estimator_.predict_log_proba(self.transform(X))

    def _more_tags(self):
        return {'poor_score': True}


class RFECV(RFE, MetaEstimatorMixin):
    """Feature ranking with recursive feature elimination and cross-validated
    selection of the best number of features.

    See glossary entry for :term:`cross-validation estimator`.

    Read more in the :ref:`User Guide <rfe>`.

    Parameters
    ----------
    estimator : object
        A supervised learning estimator with a ``fit`` method that provides
        information about feature importance either through a ``coef_``
        attribute or through a ``feature_importances_`` attribute.

    step : int or float, optional (default=1)
        If greater than or equal to 1, then ``step`` corresponds to the
        (integer) number of features to remove at each iteration.
        If within (0.0, 1.0), then ``step`` corresponds to the percentage
        (rounded down) of features to remove at each iteration.
        Note that the last iteration may remove fewer than ``step`` features in
        order to reach ``min_features_to_select``.

    step_change_at : int or float or None, optional (default=None)
        If specified ``step`` int or float equates to > 1 feature, then the
        threshold number of remaining features when to change to
        ``changed_step``. If ``step_change_at`` within (0.0, 1.0), then
        threshold number of remaining features calculated as fraction of
        original number of features.

    changed_step : int or float, optional (default=1)
        Step size to change to starting at ``step_change_at`` if specified. If
        greater than or equal to 1, then ``changed_step`` corresponds to the
        (integer) number of features to remove at each iteration. If within
        (0.0, 1.0), then ``changed_step`` corresponds to the percentage
        (rounded down) of features to remove at each iteration.

    reducing_step : boolean, optional (default=False)
        If true and ``step`` or ``changed_step`` is a float, the number of
        features removed is calculated as a fraction of the remaining features
        in that iteration. If false, the number of features removed is constant
        (a fraction of the original number of features) across iterations.

    min_features_to_select : int, (default=1)
        The minimum number of features to be selected. This number of features
        will always be scored, even if the difference between the original
        feature count and ``min_features_to_select`` isn't divisible by
        ``step``.

    cv : int, cross-validation generator or an iterable, optional
        Determines the cross-validation splitting strategy.
        Possible inputs for cv are:

        - None, to use the default 3-fold cross-validation,
        - integer, to specify the number of folds.
        - :term:`CV splitter`,
        - An iterable yielding (train, test) splits as arrays of indices.

        For integer/None inputs, if ``y`` is binary or multiclass,
        :class:`sklearn.model_selection.StratifiedKFold` is used. If the
        estimator is a classifier or if ``y`` is neither binary nor multiclass,
        :class:`sklearn.model_selection.KFold` is used.

        Refer :ref:`User Guide <cross_validation>` for the various
        cross-validation strategies that can be used here.

        .. versionchanged:: 0.20
            ``cv`` default value of None will change from 3-fold to 5-fold
            in v0.22.

    scoring : string, callable or None, optional, (default=None)
        A string (see model evaluation documentation) or
        a scorer callable object / function with signature
        ``scorer(estimator, X, y)``.

    verbose : int, (default=0)
        Controls verbosity of output.

    n_jobs : int or None, optional (default=None)
        Number of cores to run in parallel while fitting across folds.
        ``None`` means 1 unless in a :obj:`joblib.parallel_backend` context.
        ``-1`` means using all processors. See :term:`Glossary <n_jobs>`
        for more details.

    Attributes
    ----------
    n_features_ : int
        The number of selected features with cross-validation.

    support_ : array of shape [n_features]
        The mask of selected features.

    ranking_ : array of shape [n_features]
        The feature ranking, such that `ranking_[i]`
        corresponds to the ranking
        position of the i-th feature.
        Selected (i.e., estimated best)
        features are assigned rank 1.

    grid_scores_ : array of shape [n_subsets_of_features]
        The cross-validation scores such that
        ``grid_scores_[i]`` corresponds to
        the CV score of the i-th subset of features.

    estimator_ : object
        The external estimator fit on the reduced dataset.

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
    >>> selector.support_ # doctest: +NORMALIZE_WHITESPACE
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
    def __init__(self, estimator, step=1, step_change_at=None, changed_step=1,
                 reducing_step=False, min_features_to_select=1, cv='warn',
                 scoring=None, verbose=0, n_jobs=None):
        self.estimator = estimator
        self.step = step
        self.step_change_at = step_change_at
        self.changed_step = changed_step
        self.reducing_step = reducing_step
        self.cv = cv
        self.scoring = scoring
        self.verbose = verbose
        self.n_jobs = n_jobs
        self.min_features_to_select = min_features_to_select

    def fit(self, X, y, groups=None):
        """Fit the RFE model and automatically tune the number of selected
           features.

        Parameters
        ----------
        X : {array-like, sparse matrix}, shape = [n_samples, n_features]
            Training vector, where `n_samples` is the number of samples and
            `n_features` is the total number of features.

        y : array-like, shape = [n_samples]
            Target values (integers for classification, real numbers for
            regression).

        groups : array-like, shape = [n_samples], optional
            Group labels for the samples used while splitting the dataset into
            train/test set.
        """
        X, y = check_X_y(X, y, "csr", ensure_min_features=2)

        # Initialization
        cv = check_cv(self.cv, y, is_classifier(self.estimator))
        scorer = check_scoring(self.estimator, scoring=self.scoring)

        # Build an RFE object, which will evaluate and score each possible
        # feature count, down to self.min_features_to_select
        rfe = RFE(estimator=self.estimator,
                  n_features_to_select=self.min_features_to_select,
                  step=self.step, step_change_at=self.step_change_at,
                  changed_step=self.changed_step,
                  reducing_step=self.reducing_step, verbose=self.verbose)

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

        cv_results = parallel(
            func(rfe, self.estimator, X, y, train, test, scorer)
            for train, test in cv.split(X, y, groups))

        # Reverse scores and num feature steps to select argmax score with
        # lowest num features in case of a score tie
        scores = np.sum([s for s, _ in cv_results], axis=0)[::-1]
        n_remaining_feature_steps = cv_results[0][1][::-1]
        n_features_to_select = n_remaining_feature_steps[np.argmax(scores)]

        # Re-execute an elimination with best_k over the whole set
        rfe = RFE(estimator=self.estimator,
                  n_features_to_select=n_features_to_select, step=self.step,
                  step_change_at=self.step_change_at,
                  changed_step=self.changed_step,
                  reducing_step=self.reducing_step, verbose=self.verbose)

        rfe.fit(X, y)

        # Set final attributes
        self.support_ = rfe.support_
        self.n_features_ = rfe.n_features_
        self.ranking_ = rfe.ranking_
        self.estimator_ = clone(self.estimator)
        self.estimator_.fit(self.transform(X), y)

        # Fixing a normalization error, n is equal to get_n_splits(X, y) - 1
        # here, the scores are normalized by get_n_splits(X, y)
        self.grid_scores_ = scores / cv.get_n_splits(X, y, groups)
        return self
