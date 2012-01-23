"""Forest of trees-based ensemble methods

Those methods include random forests and extremely randomized trees.

The module structure is the following:

- The ``BaseForest`` base class implements a common ``fit`` method for all
  the estimators in the module. The ``fit`` method of the base ``Forest``
  class calls the ``fit`` method of each sub-estimator on random samples
  (with replacement, a.k.a. bootstrap) of the training set.

  The init of the sub-estimator is further delegated to the
  ``BaseEnsemble`` constructor.

- The ``ForestClassifier`` and ``ForestRegressor`` base classes further
  implement the prediction logic by computing an average of the predicted
  outcomes of the sub-estimators.

- The ``RandomForestClassifier`` and ``RandomForestRegressor`` derived
  classes provide the user with concrete implementations of
  the forest ensemble method using classical, deterministic
  ``DecisionTreeClassifier`` and ``DecisionTreeRegressor`` as
  sub-estimator implementations.

- The ``ExtraTreesClassifier`` and ``ExtraTreesRegressor`` derived
  classes provide the user with concrete implementations of the
  forest ensemble method using the extremly randomized trees
  ``ExtraTreeClassifier`` and ``ExtraTreeRegressor`` as
  sub-estimator implementations.

"""

# Authors: Gilles Louppe, Brian Holt
# License: BSD 3

import itertools
import numpy as np

from ..feature_selection.selector_mixin import SelectorMixin
from ..tree import DecisionTreeClassifier, DecisionTreeRegressor, \
                   ExtraTreeClassifier, ExtraTreeRegressor

from .bagging import BaggedClassifier, BaggedRegressor

__all__ = ["RandomForestClassifier",
           "RandomForestRegressor",
           "ExtraTreesClassifier",
           "ExtraTreesRegressor"]


def _pre_fit(forest, X, y):
    if forest.bootstrap:
        sample_mask = None
        X_argsorted = None

    else:
        sample_mask = np.ones((X.shape[0],), dtype=np.bool)
        X_argsorted = np.asfortranarray(
            np.argsort(X.T, axis=1).astype(np.int32).T)

    return {"X_argsorted": X_argsorted, "sample_mask": sample_mask}

def _post_fit(forest, X, y):
    if forest.compute_importances:
        forest.feature_importances_ = \
            sum(tree.feature_importances_ for tree in forest.estimators_) \
            / forest.n_estimators


class ForestClassifier(BaggedClassifier, SelectorMixin):
    """Base class for forest of trees-based classifiers.

    Warning: This class should not be used directly. Use derived classes
    instead.
    """
    def __init__(self, base_estimator,
                       n_estimators,
                       estimator_params,
                       bootstrap,
                       oob_score,
                       n_jobs,
                       random_state,
                       criterion,
                       max_depth,
                       min_split,
                       min_density,
                       max_features,
                       compute_importances):
        super(ForestClassifier, self).__init__(
            base_estimator=base_estimator,
            n_estimators=n_estimators,
            estimator_params=estimator_params,
            bootstrap=bootstrap,
            oob_score=oob_score,
            n_jobs=n_jobs,
            random_state=random_state)

        self.criterion = criterion
        self.max_depth = max_depth
        self.min_split = min_split
        self.min_density = min_density
        self.max_features = max_features
        self.compute_importances = compute_importances

        self.feature_importances_ = None

    def _pre_fit(self, X, y):
        return _pre_fit(self, X, y)

    def _post_fit(self, X, y):
        return _post_fit(self, X, y)


class ForestRegressor(BaggedRegressor, SelectorMixin):
    """Base class for forest of trees-based regressors.

    Warning: This class should not be used directly. Use derived classes
    instead.
    """
    def __init__(self, base_estimator,
                       n_estimators,
                       estimator_params,
                       bootstrap,
                       oob_score,
                       n_jobs,
                       random_state,
                       criterion,
                       max_depth,
                       min_split,
                       min_density,
                       max_features,
                       compute_importances):
        super(ForestRegressor, self).__init__(
            base_estimator=base_estimator,
            n_estimators=n_estimators,
            estimator_params=estimator_params,
            bootstrap=bootstrap,
            oob_score=oob_score,
            n_jobs=n_jobs,
            random_state=random_state)

        self.criterion = criterion
        self.max_depth = max_depth
        self.min_split = min_split
        self.min_density = min_density
        self.max_features = max_features
        self.compute_importances = compute_importances

        self.feature_importances_ = None

    def _pre_fit(self, X, y):
        return _pre_fit(self, X, y)

    def _post_fit(self, X, y):
        return _post_fit(self, X, y)


class RandomForestClassifier(ForestClassifier):
    """A random forest classifier.

    A random forest is a meta estimator that fits a number of classifical
    decision trees on various sub-samples of the dataset and use averaging
    to improve the predictive accuracy and control over-fitting.

    Parameters
    ----------
    n_estimators : integer, optional (default=10)
        The number of trees in the forest.

    criterion : string, optional (default="gini")
        The function to measure the quality of a split. Supported criteria are
        "gini" for the Gini impurity and "entropy" for the information gain.
        Note: this parameter is tree-specific.

    max_depth : integer or None, optional (default=None)
        The maximum depth of the tree. If None, then nodes are expanded until
        all leaves are pure or until all leaves contain less than min_split
        samples.
        Note: this parameter is tree-specific.

    min_split : integer, optional (default=1)
        The minimum number of samples required to split an internal node.
        Note: this parameter is tree-specific.

    min_density : float, optional (default=0.1)
        This parameter controls a trade-off in an optimization heuristic. It
        controls the minimum density of the `sample_mask` (i.e. the
        fraction of samples in the mask). If the density falls below this
        threshold the mask is recomputed and the input data is packed
        which results in data copying.  If `min_density` equals to one,
        the partitions are always represented as copies of the original
        data. Otherwise, partitions are represented as bit masks (aka
        sample masks).
        Note: this parameter is tree-specific.

    max_features : int, string or None, optional (default="auto")
        The number of features to consider when looking for the best split:
          - If "auto", then `max_features=sqrt(n_features)` on
            classification tasks and `max_features=n_features` on regression
            problems.
          - If "sqrt", then `max_features=sqrt(n_features)`.
          - If "log2", then `max_features=log2(n_features)`.
          - If None, then `max_features=n_features`.

        Note: this parameter is tree-specific.

    bootstrap : boolean, optional (default=True)
        Whether bootstrap samples are used when building trees.

    oob_score : bool
        Whether to use out-of-bag samples to estimate
        the generalization error.

    compute_importances : boolean, optional (default=True)
        Whether feature importances are computed and stored into the
        ``feature_importances_`` attribute when calling fit.

    n_jobs : integer, optional (default=1)
        The number of jobs to run in parallel. If -1, then the number of jobs
        is set to the number of cores.

    random_state : int, RandomState instance or None, optional (default=None)
        If int, random_state is the seed used by the random number generator;
        If RandomState instance, random_state is the random number generator;
        If None, the random number generator is the RandomState instance used
        by `np.random`.

    Attributes
    ----------
    `feature_importances_` : array, shape = [n_features]
        The feature importances (the higher, the more important the feature).

    `oob_score_` : float
        Score of the training dataset obtained using an out-of-bag estimate.

    `oob_decision_function_` : array, shape = [n_samples, n_classes]
        Decision function computed with out-of-bag estimate on the training
        set.

    Notes
    -----
    **References**:

    .. [1] L. Breiman, "Random Forests", Machine Learning, 45(1), 5-32, 2001.

    See also
    --------
    DecisionTreeClassifier, ExtraTreesClassifier
    """
    def __init__(self, n_estimators=10,
                       criterion="gini",
                       max_depth=None,
                       min_split=1,
                       min_density=0.1,
                       max_features="auto",
                       bootstrap=True,
                       oob_score=False,
                       compute_importances=False,
                       n_jobs=1,
                       random_state=None):
        super(RandomForestClassifier, self).__init__(
            base_estimator=DecisionTreeClassifier(),
            n_estimators=n_estimators,
            estimator_params=("criterion", "max_depth", "min_split",
                              "min_density", "max_features",
                              "compute_importances", "random_state"),
            bootstrap=bootstrap,
            oob_score=oob_score,
            n_jobs=n_jobs,
            random_state=random_state,
            criterion=criterion,
            max_depth=max_depth,
            min_split=min_split,
            min_density=min_density,
            max_features=max_features,
            compute_importances=compute_importances)


class RandomForestRegressor(ForestRegressor):
    """A random forest regressor.

    A random forest is a meta estimator that fits a number of classifical
    decision trees on various sub-samples of the dataset and use averaging
    to improve the predictive accuracy and control over-fitting.

    Parameters
    ----------
    n_estimators : integer, optional (default=10)
        The number of trees in the forest.

    criterion : string, optional (default="mse")
        The function to measure the quality of a split. The only supported
        criterion is "mse" for the mean squared error.
        Note: this parameter is tree-specific.

    max_depth : integer or None, optional (default=None)
        The maximum depth of the tree. If None, then nodes are expanded until
        all leaves are pure or until all leaves contain less than min_split
        samples.
        Note: this parameter is tree-specific.

    min_split : integer, optional (default=1)
        The minimum number of samples required to split an internal node.
        Note: this parameter is tree-specific.

    min_density : float, optional (default=0.1)
        This parameter controls a trade-off in an optimization heuristic. It
        controls the minimum density of the `sample_mask` (i.e. the
        fraction of samples in the mask). If the density falls below this
        threshold the mask is recomputed and the input data is packed
        which results in data copying.  If `min_density` equals to one,
        the partitions are always represented as copies of the original
        data. Otherwise, partitions are represented as bit masks (aka
        sample masks).
        Note: this parameter is tree-specific.

    max_features : int, string or None, optional (default="auto")
        The number of features to consider when looking for the best split:
          - If "auto", then `max_features=sqrt(n_features)` on
            classification tasks and `max_features=n_features`
            on regression problems.
          - If "sqrt", then `max_features=sqrt(n_features)`.
          - If "log2", then `max_features=log2(n_features)`.
          - If None, then `max_features=n_features`.

        Note: this parameter is tree-specific.

    bootstrap : boolean, optional (default=True)
        Whether bootstrap samples are used when building trees.

    oob_score : bool
        whether to use out-of-bag samples to estimate
        the generalization error.

    compute_importances : boolean, optional (default=True)
        Whether feature importances are computed and stored into the
        ``feature_importances_`` attribute when calling fit.

    n_jobs : integer, optional (default=1)
        The number of jobs to run in parallel. If -1, then the number of jobs
        is set to the number of cores.

    random_state : int, RandomState instance or None, optional (default=None)
        If int, random_state is the seed used by the random number generator;
        If RandomState instance, random_state is the random number generator;
        If None, the random number generator is the RandomState instance used
        by `np.random`.

    Attributes
    ----------
    `feature_importances_` : array of shape = [n_features]
        The feature mportances (the higher, the more important the feature).

    `oob_score_` : float
        Score of the training dataset obtained using an out-of-bag estimate.

    `oob_prediction_` : array, shape = [n_samples]
        Prediction computed with out-of-bag estimate on the training set.

    Notes
    -----
    **References**:

    .. [1] L. Breiman, "Random Forests", Machine Learning, 45(1), 5-32, 2001.

    See also
    --------
    DecisionTreeRegressor, ExtraTreesRegressor
    """
    def __init__(self, n_estimators=10,
                       criterion="mse",
                       max_depth=None,
                       min_split=1,
                       min_density=0.1,
                       max_features="auto",
                       bootstrap=True,
                       oob_score=False,
                       compute_importances=False,
                       n_jobs=1,
                       random_state=None):
        super(RandomForestRegressor, self).__init__(
            base_estimator=DecisionTreeRegressor(),
            n_estimators=n_estimators,
            estimator_params=("criterion", "max_depth", "min_split",
                              "min_density", "max_features",
                              "compute_importances", "random_state"),
            bootstrap=bootstrap,
            oob_score=oob_score,
            n_jobs=n_jobs,
            random_state=random_state,
            criterion=criterion,
            max_depth=max_depth,
            min_split=min_split,
            min_density=min_density,
            max_features=max_features,
            compute_importances=compute_importances)


class ExtraTreesClassifier(ForestClassifier):
    """An extra-trees classifier.

    This class implements a meta estimator that fits a number of
    randomized decision trees (a.k.a. extra-trees) on various sub-samples
    of the dataset and use averaging to improve the predictive accuracy
    and control over-fitting.

    Parameters
    ----------
    n_estimators : integer, optional (default=10)
        The number of trees in the forest.

    criterion : string, optional (default="gini")
        The function to measure the quality of a split. Supported criteria are
        "gini" for the Gini impurity and "entropy" for the information gain.
        Note: this parameter is tree-specific.

    max_depth : integer or None, optional (default=None)
        The maximum depth of the tree. If None, then nodes are expanded until
        all leaves are pure or until all leaves contain less than min_split
        samples.
        Note: this parameter is tree-specific.

    min_split : integer, optional (default=1)
        The minimum number of samples required to split an internal node.
        Note: this parameter is tree-specific.

    min_density : float, optional (default=0.1)
        This parameter controls a trade-off in an optimization heuristic. It
        controls the minimum density of the `sample_mask` (i.e. the
        fraction of samples in the mask). If the density falls below this
        threshold the mask is recomputed and the input data is packed
        which results in data copying.  If `min_density` equals to one,
        the partitions are always represented as copies of the original
        data. Otherwise, partitions are represented as bit masks (aka
        sample masks).
        Note: this parameter is tree-specific.

    max_features : int, string or None, optional (default="auto")
        The number of features to consider when looking for the best split.
          - If "auto", then `max_features=sqrt(n_features)` on
            classification tasks and `max_features=n_features`
            on regression problems.
          - If "sqrt", then `max_features=sqrt(n_features)`.
          - If "log2", then `max_features=log2(n_features)`.
          - If None, then `max_features=n_features`.

        Note: this parameter is tree-specific.

    bootstrap : boolean, optional (default=False)
        Whether bootstrap samples are used when building trees.

    oob_score : bool
        Whether to use out-of-bag samples to estimate
        the generalization error.

    compute_importances : boolean, optional (default=True)
        Whether feature importances are computed and stored into the
        ``feature_importances_`` attribute when calling fit.

    n_jobs : integer, optional (default=1)
        The number of jobs to run in parallel. If -1, then the number of jobs
        is set to the number of cores.

    random_state : int, RandomState instance or None, optional (default=None)
        If int, random_state is the seed used by the random number generator;
        If RandomState instance, random_state is the random number generator;
        If None, the random number generator is the RandomState instance used
        by `np.random`.

    Attributes
    ----------
    `feature_importances_` : array of shape = [n_features]
        The feature mportances (the higher, the more important the feature).

    `oob_score_` : float
        Score of the training dataset obtained using an out-of-bag estimate.

    `oob_decision_function_` : array, shape = [n_samples, n_classes]
        Decision function computed with out-of-bag estimate on the training
        set.

    Notes
    -----
    **References**:

    .. [1] P. Geurts, D. Ernst., and L. Wehenkel, "Extremely randomized trees",
           Machine Learning, 63(1), 3-42, 2006.

    See also
    --------
    sklearn.tree.ExtraTreeClassifier : Base classifier for this ensemble.
    RandomForestClassifier : Ensemble Classifier based on trees with optimal
        splits.
    """
    def __init__(self, n_estimators=10,
                       criterion="gini",
                       max_depth=None,
                       min_split=1,
                       min_density=0.1,
                       max_features="auto",
                       bootstrap=False,
                       oob_score=False,
                       compute_importances=False,
                       n_jobs=1,
                       random_state=None):
        super(ExtraTreesClassifier, self).__init__(
            base_estimator=ExtraTreeClassifier(),
            n_estimators=n_estimators,
            estimator_params=("criterion", "max_depth", "min_split",
                              "min_density", "max_features",
                              "compute_importances", "random_state"),
            bootstrap=bootstrap,
            oob_score=oob_score,
            n_jobs=n_jobs,
            random_state=random_state,
            criterion=criterion,
            max_depth=max_depth,
            min_split=min_split,
            min_density=min_density,
            max_features=max_features,
            compute_importances=compute_importances)


class ExtraTreesRegressor(ForestRegressor):
    """An extra-trees regressor.

    This class implements a meta estimator that fits a number of
    randomized decision trees (a.k.a. extra-trees) on various sub-samples
    of the dataset and use averaging to improve the predictive accuracy
    and control over-fitting.

    Parameters
    ----------
    n_estimators : integer, optional (default=10)
        The number of trees in the forest.

    criterion : string, optional (default="mse")
        The function to measure the quality of a split. The only supported
        criterion is "mse" for the mean squared error.
        Note: this parameter is tree-specific.

    max_depth : integer or None, optional (default=None)
        The maximum depth of the tree. If None, then nodes are expanded until
        all leaves are pure or until all leaves contain less than min_split
        samples.
        Note: this parameter is tree-specific.

    min_split : integer, optional (default=1)
        The minimum number of samples required to split an internal node.
        Note: this parameter is tree-specific.

    min_density : float, optional (default=0.1)
        This parameter controls a trade-off in an optimization heuristic. It
        controls the minimum density of the `sample_mask` (i.e. the
        fraction of samples in the mask). If the density falls below this
        threshold the mask is recomputed and the input data is packed
        which results in data copying.  If `min_density` equals to one,
        the partitions are always represented as copies of the original
        data. Otherwise, partitions are represented as bit masks (aka
        sample masks).
        Note: this parameter is tree-specific.

    max_features : int, string or None, optional (default="auto")
        The number of features to consider when looking for the best split:
          - If "auto", then `max_features=sqrt(n_features)` on
            classification tasks and `max_features=n_features`
            on regression problems.
          - If "sqrt", then `max_features=sqrt(n_features)`.
          - If "log2", then `max_features=log2(n_features)`.
          - If None, then `max_features=n_features`.

        Note: this parameter is tree-specific.

    bootstrap : boolean, optional (default=False)
        Whether bootstrap samples are used when building trees.
        Note: this parameter is tree-specific.

    oob_score : bool
        Whether to use out-of-bag samples to estimate
        the generalization error.

    compute_importances : boolean, optional (default=True)
        Whether feature importances are computed and stored into the
        ``feature_importances_`` attribute when calling fit.

    n_jobs : integer, optional (default=1)
        The number of jobs to run in parallel. If -1, then the number of jobs
        is set to the number of cores.

    random_state : int, RandomState instance or None, optional (default=None)
        If int, random_state is the seed used by the random number generator;
        If RandomState instance, random_state is the random number generator;
        If None, the random number generator is the RandomState instance used
        by `np.random`.

    Attributes
    ----------
    `feature_importances_` : array of shape = [n_features]
        The feature mportances (the higher, the more important the feature).

    `oob_score_` : float
        Score of the training dataset obtained using an out-of-bag estimate.

    `oob_prediction_` : array, shape = [n_samples]
        Prediction computed with out-of-bag estimate on the training set.

    Notes
    -----
    **References**:

    .. [1] P. Geurts, D. Ernst., and L. Wehenkel, "Extremely randomized trees",
           Machine Learning, 63(1), 3-42, 2006.

    See also
    --------
    sklearn.tree.ExtraTreeRegressor: Base estimator for this ensemble.
    RandomForestRegressor: Ensemble regressor using trees with optimal splits.
    """
    def __init__(self, n_estimators=10,
                       criterion="mse",
                       max_depth=None,
                       min_split=1,
                       min_density=0.1,
                       max_features="auto",
                       bootstrap=False,
                       oob_score=False,
                       compute_importances=False,
                       n_jobs=1,
                       random_state=None):
        super(ExtraTreesRegressor, self).__init__(
            base_estimator=ExtraTreeRegressor(),
            n_estimators=n_estimators,
            estimator_params=("criterion", "max_depth", "min_split",
                              "min_density", "max_features",
                              "compute_importances", "random_state"),
            bootstrap=bootstrap,
            oob_score=oob_score,
            n_jobs=n_jobs,
            random_state=random_state,
            criterion=criterion,
            max_depth=max_depth,
            min_split=min_split,
            min_density=min_density,
            max_features=max_features,
            compute_importances=compute_importances)
