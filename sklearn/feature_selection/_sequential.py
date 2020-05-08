"""
Sequential feature selection
"""

import numpy as np

from .base import SelectorMixin
from ..base import BaseEstimator
from ..base import MetaEstimatorMixin
from ..base import clone
from ..utils.validation import check_is_fitted
from ..utils import _safe_indexing
from ..model_selection import cross_val_score


class SequentialFeatureSelector(SelectorMixin, MetaEstimatorMixin,
                                BaseEstimator):
    """Transformer that performs Sequential Feature Selection.

    This Sequential Feature Selector adds (forward selection) or
    removes (backward selection) features to form a feature subset in a
    greedy fashion. At each stage, this estimator chooses the best feature to
    add or remove based on the cross-validation score of an estimator.

    Read more in the :ref:`User Guide <sequential_feature_selection>`.

    .. versionadded:: 0.24

    Parameters
    ----------
    estimator : estimator instance
        An unfitted estimator

    n_features_to_select : int, default=None
        The number of features to select. If None, half of the features
        are selected.

    forward : bool, default=True
        Whether to perform forward selection or backward selection

    scoring : str, callable, list/tuple or dict, default=None
        A single str (see :ref:`scoring_parameter`) or a callable
        (see :ref:`scoring`) to evaluate the predictions on the test set.

        NOTE that when using custom scorers, each scorer should return a single
        value. Metric functions returning a list/array of values can be wrapped
        into multiple scorers that return one value each.

        If None, the estimator's score method is used.

    cv : int, cross-validation generator or an iterable, default=None
        Determines the cross-validation splitting strategy.
        Possible inputs for cv are:

        - None, to use the default 5-fold cross validation,
        - integer, to specify the number of folds in a `(Stratified)KFold`,
        - :term:`CV splitter`,
        - An iterable yielding (train, test) splits as arrays of indices.

        For integer/None inputs, if the estimator is a classifier and ``y`` is
        either binary or multiclass, :class:`StratifiedKFold` is used. In all
        other cases, :class:`KFold` is used.

        Refer :ref:`User Guide <cross_validation>` for the various
        cross-validation strategies that can be used here.

    n_jobs : int, default=None
        Number of jobs to run in parallel.
        ``None`` means 1 unless in a :obj:`joblib.parallel_backend` context.
        ``-1`` means using all processors. See :term:`Glossary <n_jobs>`
        for more details.

    pre_dispatch : int, or str, default=None
        Controls the number of jobs that get dispatched during parallel
        execution. Reducing this number can be useful to avoid an
        explosion of memory consumption when more jobs get dispatched
        than CPUs can process. This parameter can be:

            - None, in which case all the jobs are immediately
              created and spawned. Use this for lightweight and
              fast-running jobs, to avoid delays due to on-demand
              spawning of the jobs

            - An int, giving the exact number of total jobs that are
              spawned

            - A str, giving an expression as a function of n_jobs,
              as in '2*n_jobs'

    Attributes
    ----------
    n_features_to_select_ : int
        The number of features that were selected. It corresponds to
        `n_features_to_select` unless the parameter was None.

    support_ : ndarray of bool of shape (n_features,)
        The mask of selected features.

    Examples
    --------
    >>> from sklearn.feature_selection import SequentialFeatureSelector
    >>> from sklearn.neighbors import KNeighborsClassifier
    >>> from sklearn.datasets import load_iris
    >>> X, y = load_iris(return_X_y=True)
    >>> knn = KNeighborsClassifier(n_neighbors=3)
    >>> sfs = SequentialFeatureSelector(knn, n_features_to_select=3)
    >>> sfs.fit(X, y)
    SequentialFeatureSelector(estimator=KNeighborsClassifier(n_neighbors=3),
                              n_features_to_select=3)
    >>> sfs.get_support()
    array([ True, False,  True,  True])
    >>> sfs.transform(X).shape
    (150, 3)

    See also
    --------
    RFE : Recursive feature elimination selection of the best number of
          features
    SelectFromModel : Feature selection based on importance weights
    """
    def __init__(self, estimator, n_features_to_select=None, forward=True,
                 scoring=None, cv=5, n_jobs=None, pre_dispatch='2*n_jobs'):

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
        X : array-like of shape (n_samples, n_features)
            Training vectors.
        y : array-like of shape (n_samples,)
            Target values.

        Returns
        -------
        self : object
        """

        tags = self._get_tags()
        X, y = self._validate_data(
            X, y, accept_sparse="csc",
            ensure_min_features=2,
            force_all_finite=not tags.get('allow_nan', True),
            multi_output=True
        )

        if self.n_features_to_select is None:
            self.n_features_to_select_ = X.shape[1] // 2
        else:
            self.n_features_to_select_ = self.n_features_to_select

        if not 1 <= self.n_features_to_select_ < X.shape[1]:
            raise ValueError(
                "n_features_to_select must be in [1, n_features - 1] = "
                f"[1, {X.shape[1] - 1}]. Got {self.n_features_to_select_}."
            )

        cloned_estimator = clone(self.estimator)

        # the current mask corresponds to the set of features:
        # - that we have already selected when doing forward selection
        # - that we have already excluded when doing backward selection
        current_mask = set()
        n_iterations = (self.n_features_to_select_ if self.forward
                        else X.shape[1] - self.n_features_to_select_)
        for _ in range(n_iterations):
            new_feature_idx = self._get_best_new_feature(cloned_estimator, X,
                                                         y, current_mask)
            current_mask.add(new_feature_idx)

        # transform the mask into a proper boolean mask of selected features
        self.support_ = np.full(X.shape[1], fill_value=not self.forward,
                                dtype=bool)
        self.support_[list(current_mask)] = self.forward

        return self

    def _get_best_new_feature(self, estimator, X, y, current_mask):
        # Return the best new feature to add to the current_mask, i.e. return
        # the best new feature to add (resp. remove) when doing forward
        # selection (resp. backward selection)
        all_features = set(range(X.shape[1]))
        candidate_feature_indices = all_features - current_mask
        scores = {}
        for feature_idx in candidate_feature_indices:
            new_mask = current_mask | {feature_idx}
            if not self.forward:
                # transform new_mask into its complement, i.e. change its
                # semantic from "features to remove" into "features to keep"
                # because _safe_indexing only understands the latter
                # TODO: remove when _safe_indexing supports "complement"
                new_mask = all_features - new_mask
            X_new = _safe_indexing(X, list(new_mask), axis=1)
            scores[feature_idx] = cross_val_score(
                estimator, X_new, y, cv=self.cv, scoring=self.scoring,
                n_jobs=self.n_jobs, pre_dispatch=self.pre_dispatch).mean()
        return max(scores, key=lambda feature_idx: scores[feature_idx])

    def _get_support_mask(self):
        check_is_fitted(self)
        return self.support_

    def _more_tags(self):
        estimator_tags = self.estimator._get_tags()
        return {'allow_nan': estimator_tags.get('allow_nan', True),
                'requires_y': True,
                }
