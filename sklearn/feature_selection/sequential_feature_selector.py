"""
Sequential feature selection

"""

# Author: Sebastian Raschka <se.raschka@gmail.com>
#
# License: BSD 3 clause

import numpy as np
from itertools import combinations
from .base import SelectorMixin
from ..base import BaseEstimator
from ..base import MetaEstimatorMixin
from ..base import clone
from ..utils.validation import check_is_fitted
from ..model_selection import cross_val_score
from ..metrics.scorer import check_scoring


class SequentialFeatureSelector(BaseEstimator, MetaEstimatorMixin,
                                SelectorMixin):
    """Feature selector that selects features via greedy search.

    This Sequential Feature Selector adds (forward selection) or
    removes (backward selection) the features (X) to form a feature subset
    in a greedy fashion that optimizes the extrinsic performance metric
    of a Regressor or Classifier on the desired ouputs (y).

    Read more in the :ref:`User Guide <sequential_feature_selection>`.

    Parameters
    ----------
    estimator : scikit-learn Classifier or Regressor
        Invoking the ``fit`` method on the `SequentialFeatureSelector`
        will fit a clone of the original estimator that
        will be stored in the class attribute `self.estimators_`.

    n_features_to_select : int or tuple (default=1)
        An integer arguments specifies the number of features to select,
        where n_features_to_select < the full feature set.
        Optionally, a tuple containing a min and max value can be provided
        so that the feature selector will return a feature subset between
        with min <= n_features <= max that scored highest in the evaluation.
        For example, the tuple (1, 4) will return any combination from
        1 up to 4 features instead of a fixed number
        of features `n_features_to_select`.

    scoring : str, callable, or None (default=None)
        A string (see model evaluation documentation) or a scorer
        callable object / function with signature `scorer(estimator, X, y)`.

    forward : bool (default=True)
        Performs forward selection if True and backward selection, otherwise.

    cv : int or cross-validation generator, or an iterable (default=5)
        Determines the cross-validation splitting strategy for
        feature selection. Possible inputs for cv are:
        - 0 or None, don't use cross validation
        - integer > 1, to specify the number of folds in a (Stratified)KFold
        - An object to be used as a cross-validation generator.
        - An iterable yielding train, test splits.
        For integer/None inputs, if the estimator is a classifier
        and `y` is either binary or multiclass, `StratifiedKFold` is used.
        In all other cases, `KFold` is used.

    n_jobs : int (default=1)
        The number of CPUs to use for cross validation.

    pre_dispatch : int, or string (default: '2*n_jobs')
        Controls the number of jobs that get dispatched
        during parallel execution in cross_val_score.
        Reducing this number can be useful to avoid an explosion of
        memory consumption when more jobs get dispatched than CPUs can process.
        This parameter can be:
        - None, in which case all the jobs are immediately created and spawned.
          Use this for lightweight and fast-running jobs,
          to avoid delays due to on-demand spawning of the jobs
        - An int, giving the exact number of total jobs that are spawned
        - A string, giving an expression as a function
            of `n_jobs`, as in `2*n_jobs`

    Attributes
    ----------
    estimator_ : scikit-learn Classifier or Regressor
        Fitted clone of `estimator`.

    feature_subset_idx_ : array-like, shape = [n_predictions]
        Feature Indices of the selected feature subsets.

    score_ : float
        Cross validation average score of the selected subset.

    subsets_ : dict
        A dictionary of selected feature subsets during the
        sequential selection, where the dictionary keys are
        the lengths k of these feature subsets. The dictionary
        values are dictionaries themselves with the following
        keys: 'feature_subset_idx' (tuple of indices of the feature subset)
              'cv_scores' (list individual cross-validation scores)
              'avg_score' (average cross-validation score)

    Examples
    --------
    The following example shows how to use the sequential feature selector
    with default settings to select a feature subset, consisting of
    1 to 3 features, from iris. The selection criteria for this
    feature subset is the average cross-validation performance
    (cv=5 by default) of the `estimator` (here: KNN)
    during the greedy forward selection search.

        >>> from sklearn.feature_selection import SequentialFeatureSelector
        >>> from sklearn.neighbors import KNeighborsClassifier
        >>> from sklearn.datasets import load_iris
        >>> iris = load_iris()
        >>> X, y = iris.data, iris.target
        >>> knn = KNeighborsClassifier(n_neighbors=3)
        >>> sfs = SequentialFeatureSelector(knn, n_features_to_select=(1, 3))
        >>> sfs = sfs.fit(X, y)
        >>> sfs.score_  # doctest: +ELLIPSIS
        0.9733...
        >>> sfs.feature_subset_idx_
        (0, 2, 3)
        >>> sfs.transform(X).shape
        (150, 3)

    """
    def __init__(self, estimator, n_features_to_select=1,
                 forward=True, scoring=None,
                 cv=5, n_jobs=1,
                 pre_dispatch='2*n_jobs'):

        self.estimator = estimator
        self.n_features_to_select = n_features_to_select
        self.forward = forward
        self.pre_dispatch = pre_dispatch
        self.scoring = scoring
        self.cv = cv
        self.n_jobs = n_jobs

    def fit(self, X, y):
        """Perform feature selection and learn model from training data.

        Parameters
        ----------
        X : array-like, shape = [n_samples, n_features]
            Training vectors, where n_samples is the number of samples and
            n_features is the number of features.
        y : array-like, shape = [n_samples]
            Target values.

        Returns
        -------
        self : object

        """
        if self.scoring is None:
            self.scorer = check_scoring(self.estimator)
        else:
            self.scorer = self.scoring

        if not isinstance(self.n_features_to_select, int) and\
                not isinstance(self.n_features_to_select, tuple):
            raise ValueError('n_features_to_select must be a positive integer'
                             ' or tuple')

        if isinstance(self.n_features_to_select, int) and (
                    self.n_features_to_select < 1 or
                    self.n_features_to_select > X.shape[1]):
            raise ValueError('n_features_to_select must be a positive integer'
                             ' between 1 and X.shape[1], got %s'
                             % (self.n_features_to_select, ))

        if isinstance(self.n_features_to_select, tuple):
            if len(self.n_features_to_select) != 2:
                raise ValueError('n_features_to_select tuple must consist of 2'
                                 ' elements a min and a max value.')

            if self.n_features_to_select[0] not in range(1, X.shape[1] + 1):
                raise ValueError('n_features_to_select tuple min value must'
                                 ' be in range(1, X.shape[1]+1).')

            if self.n_features_to_select[1] not in range(1, X.shape[1] + 1):
                raise ValueError('n_features_to_select tuple max value must'
                                 ' be in range(1, X.shape[1]+1).')

            if self.n_features_to_select[0] > self.n_features_to_select[1]:
                raise ValueError('The min n_features_to_select value must be'
                                 ' larger than the max'
                                 ' n_features_to_select value.')

        if isinstance(self.n_features_to_select, tuple):
            select_in_range = True
        else:
            select_in_range = False
            k_to_select = self.n_features_to_select

        self.estimator_ = clone(self.estimator)

        self._n_features = X.shape[1]
        self.subsets_ = {}
        orig_set = set(range(self._n_features))
        if self.forward:
            if select_in_range:
                k_to_select = self.n_features_to_select[1]
            k_idx = ()
            k = 0
        else:
            if select_in_range:
                k_to_select = self.n_features_to_select[0]
            k_idx = tuple(range(self._n_features))
            k = len(k_idx)
            k_score = self._calc_score(X, y, k_idx)
            self.subsets_[k] = {
                'feature_subset_idx': k_idx,
                'cv_scores': k_score,
                'avg_score': k_score.mean()
                }

        best_subset = None
        k_score = 0

        while k != k_to_select:
            prev_subset = set(k_idx)
            if self.forward:
                k_idx, k_score, cv_scores = self._inclusion(
                    orig_set=orig_set,
                    subset=prev_subset,
                    X=X,
                    y=y
                )
            else:
                k_idx, k_score, cv_scores = self._exclusion(
                    feature_set=prev_subset,
                    X=X,
                    y=y
                )

            k = len(k_idx)
            if k not in self.subsets_ or (self.subsets_[k]['avg_score'] <
                                          k_score):
                self.subsets_[k] = {
                    'feature_subset_idx': k_idx,
                    'cv_scores': cv_scores,
                    'avg_score': k_score
                }

        if select_in_range:
            max_score = float('-inf')
            for k in self.subsets_:
                if self.subsets_[k]['avg_score'] > max_score:
                    max_score = self.subsets_[k]['avg_score']
                    best_subset = k
            k_score = max_score
            k_idx = self.subsets_[best_subset]['feature_subset_idx']

        self.feature_subset_idx_ = k_idx
        self.score_ = k_score
        return self

    def _calc_score(self, X, y, indices):
        if self.cv:
            scores = cross_val_score(self.estimator_,
                                     X[:, indices], y,
                                     cv=self.cv,
                                     scoring=self.scorer,
                                     n_jobs=self.n_jobs,
                                     pre_dispatch=self.pre_dispatch)
        else:
            self.estimator_.fit(X[:, indices], y)
            scores = np.array([self.scorer(self.estimator_, X[:, indices], y)])
        return scores

    def _inclusion(self, orig_set, subset, X, y):
        all_avg_scores = []
        all_cv_scores = []
        all_subsets = []
        res = (None, None, None)
        remaining = orig_set - subset
        if remaining:
            for feature in remaining:
                new_subset = tuple(subset | {feature})
                cv_scores = self._calc_score(X, y, new_subset)
                all_avg_scores.append(cv_scores.mean())
                all_cv_scores.append(cv_scores)
                all_subsets.append(new_subset)
            best = np.argmax(all_avg_scores)
            res = (all_subsets[best],
                   all_avg_scores[best],
                   all_cv_scores[best])
        return res

    def _exclusion(self, feature_set, X, y, fixed_feature=None):
        n = len(feature_set)
        res = (None, None, None)
        if n > 1:
            all_avg_scores = []
            all_cv_scores = []
            all_subsets = []
            for p in combinations(feature_set, r=n - 1):
                if fixed_feature and fixed_feature not in set(p):
                    continue
                cv_scores = self._calc_score(X, y, p)
                all_avg_scores.append(cv_scores.mean())
                all_cv_scores.append(cv_scores)
                all_subsets.append(p)
            best = np.argmax(all_avg_scores)
            res = (all_subsets[best],
                   all_avg_scores[best],
                   all_cv_scores[best])
        return res

    @property
    def _estimator_type(self):
        return self.estimator_._estimator_type

    def _get_support_mask(self):
        check_is_fitted(self, 'feature_subset_idx_')
        mask = np.zeros((self._n_features,), dtype=np.bool)
        # list to avoid IndexError in old NumPy versions
        mask[list(self.feature_subset_idx_)] = True
        return mask
