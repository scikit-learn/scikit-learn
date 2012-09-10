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

Single and multi-output problems are both handled.

"""

# Authors: Gilles Louppe, Brian Holt
# License: BSD 3

import itertools
import numpy as np
from warnings import warn
from abc import ABCMeta, abstractmethod

from ..base import ClassifierMixin, RegressorMixin
from ..externals.joblib import Parallel, delayed, cpu_count
from ..feature_selection.selector_mixin import SelectorMixin
from ..tree import DecisionTreeClassifier, DecisionTreeRegressor, \
                   ExtraTreeClassifier, ExtraTreeRegressor
from ..tree._tree import DTYPE, DOUBLE
from ..utils import array2d, check_random_state, check_arrays
from ..metrics import r2_score

from .base import BaseEnsemble

__all__ = ["RandomForestClassifier",
           "RandomForestRegressor",
           "ExtraTreesClassifier",
           "ExtraTreesRegressor"]

MAX_INT = np.iinfo(np.int32).max


def _parallel_build_trees(n_trees, forest, X, y,
                          sample_mask, X_argsorted, seed, verbose):
    """Private function used to build a batch of trees within a job."""
    random_state = check_random_state(seed)
    trees = []

    for i in xrange(n_trees):
        if verbose > 1:
            print("building tree %d of %d" % (i + 1, n_trees))
        seed = random_state.randint(MAX_INT)

        tree = forest._make_estimator(append=False)
        tree.set_params(compute_importances=forest.compute_importances)
        tree.set_params(random_state=check_random_state(seed))

        if forest.bootstrap:
            n_samples = X.shape[0]
            indices = random_state.randint(0, n_samples, n_samples)
            tree.fit(X[indices], y[indices],
                     sample_mask=sample_mask, X_argsorted=X_argsorted)
            tree.indices_ = indices

        else:
            tree.fit(X, y,
                     sample_mask=sample_mask, X_argsorted=X_argsorted)

        trees.append(tree)

    return trees


def _parallel_predict_proba(trees, X, n_classes, n_outputs):
    """Private function used to compute a batch of predictions within a job."""
    n_samples = X.shape[0]
    p = []

    for k in xrange(n_outputs):
        p.append(np.zeros((n_samples, n_classes[k])))

    for tree in trees:
        p_tree = tree.predict_proba(X)

        if n_outputs == 1:
            p_tree = [p_tree]

        for k in xrange(n_outputs):
            if n_classes[k] == tree.n_classes_[k]:
                p[k] += p_tree[k]

            else:
                for j, c in enumerate(tree.classes_[k]):
                    p[k][:, c] += p_tree[k][:, j]

    return p


def _parallel_predict_regression(trees, X):
    """Private function used to compute a batch of predictions within a job."""
    return sum(tree.predict(X) for tree in trees)


def _partition_trees(forest):
    """Private function used to partition trees between jobs."""
    # Compute the number of jobs
    if forest.n_jobs == -1:
        n_jobs = min(cpu_count(), forest.n_estimators)

    else:
        n_jobs = min(forest.n_jobs, forest.n_estimators)

    # Partition trees between jobs
    n_trees = [int(forest.n_estimators / n_jobs)] * n_jobs

    for i in xrange(forest.n_estimators % n_jobs):
        n_trees[i] += 1

    starts = [0] * (n_jobs + 1)

    for i in xrange(1, n_jobs + 1):
        starts[i] = starts[i - 1] + n_trees[i - 1]

    return n_jobs, n_trees, starts


def _parallel_X_argsort(X):
    """Private function used to sort the features of X."""
    return np.asarray(np.argsort(X.T, axis=1).T, dtype=np.int32, order="F")


def _partition_features(forest, n_total_features):
    """Private function used to partition features between jobs."""
    # Compute the number of jobs
    if forest.n_jobs == -1:
        n_jobs = min(cpu_count(), n_total_features)

    else:
        n_jobs = min(forest.n_jobs, n_total_features)

    # Partition features between jobs
    n_features = [n_total_features / n_jobs] * n_jobs

    for i in xrange(n_total_features % n_jobs):
        n_features[i] += 1

    starts = [0] * (n_jobs + 1)

    for i in xrange(1, n_jobs + 1):
        starts[i] = starts[i - 1] + n_features[i - 1]

    return n_jobs, n_features, starts


class BaseForest(BaseEnsemble, SelectorMixin):
    """Base class for forests of trees.

    Warning: This class should not be used directly. Use derived classes
    instead.
    """
    __metaclass__ = ABCMeta

    @abstractmethod
    def __init__(self, base_estimator,
                       n_estimators=10,
                       estimator_params=[],
                       bootstrap=False,
                       compute_importances=False,
                       oob_score=False,
                       n_jobs=1,
                       random_state=None,
                       verbose=0):
        super(BaseForest, self).__init__(
            base_estimator=base_estimator,
            n_estimators=n_estimators,
            estimator_params=estimator_params)

        self.bootstrap = bootstrap
        self.compute_importances = compute_importances
        self.oob_score = oob_score
        self.n_jobs = n_jobs
        self.random_state = random_state

        self.n_features_ = None
        self.n_outputs_ = None
        self.classes_ = None
        self.n_classes_ = None
        self.feature_importances_ = None

        self.verbose = verbose

    def fit(self, X, y):
        """Build a forest of trees from the training set (X, y).

        Parameters
        ----------
        X : array-like of shape = [n_samples, n_features]
            The training input samples.

        y : array-like, shape = [n_samples] or [n_samples, n_outputs]
            The target values (integers that correspond to classes in
            classification, real numbers in regression).

        Returns
        -------
        self : object
            Returns self.
        """
        self.random_state = check_random_state(self.random_state)

        # Precompute some data
        X, y = check_arrays(X, y, sparse_format="dense")
        if getattr(X, "dtype", None) != DTYPE or \
           X.ndim != 2 or not X.flags.fortran:
            X = array2d(X, dtype=DTYPE, order="F")

        n_samples, self.n_features_ = X.shape

        if self.bootstrap:
            sample_mask = None
            X_argsorted = None

        else:
            if self.oob_score:
                raise ValueError("Out of bag estimation only available"
                        " if bootstrap=True")

            sample_mask = np.ones((n_samples,), dtype=np.bool)

            n_jobs, _, starts = _partition_features(self, self.n_features_)

            all_X_argsorted = Parallel(n_jobs=n_jobs)(
                delayed(_parallel_X_argsort)(
                    X[:, starts[i]:starts[i + 1]])
                for i in xrange(n_jobs))

            X_argsorted = np.asfortranarray(np.hstack(all_X_argsorted))

        y = np.atleast_1d(y)
        if y.ndim == 1:
            y = y[:, np.newaxis]

        self.classes_ = []
        self.n_classes_ = []
        self.n_outputs_ = y.shape[1]

        if isinstance(self.base_estimator, ClassifierMixin):
            y = np.copy(y)

            for k in xrange(self.n_outputs_):
                unique = np.unique(y[:, k])
                self.classes_.append(unique)
                self.n_classes_.append(unique.shape[0])
                y[:, k] = np.searchsorted(unique, y[:, k])

        if getattr(y, "dtype", None) != DTYPE or not y.flags.contiguous:
            y = np.ascontiguousarray(y, dtype=DOUBLE)

        # Assign chunk of trees to jobs
        n_jobs, n_trees, _ = _partition_trees(self)

        # Parallel loop
        all_trees = Parallel(n_jobs=n_jobs, verbose=self.verbose)(
            delayed(_parallel_build_trees)(
                n_trees[i],
                self,
                X,
                y,
                sample_mask,
                X_argsorted,
                self.random_state.randint(MAX_INT),
                verbose=self.verbose)
            for i in xrange(n_jobs))

        # Reduce
        self.estimators_ = [tree for tree in itertools.chain(*all_trees)]

        # Calculate out of bag predictions and score
        if self.oob_score:
            if isinstance(self, ClassifierMixin):
                self.oob_decision_function_ = []
                self.oob_score_ = 0.0

                predictions = []
                for k in xrange(self.n_outputs_):
                    predictions.append(np.zeros((n_samples,
                                                 self.n_classes_[k])))

                for estimator in self.estimators_:
                    mask = np.ones(n_samples, dtype=np.bool)
                    mask[estimator.indices_] = False

                    p_estimator = estimator.predict_proba(X[mask, :])
                    if self.n_outputs_ == 1:
                        p_estimator = [p_estimator]

                    for k in xrange(self.n_outputs_):
                        predictions[k][mask, :] += p_estimator[k]

                for k in xrange(self.n_outputs_):
                    if (predictions[k].sum(axis=1) == 0).any():
                        warn("Some inputs do not have OOB scores. "
                             "This probably means too few trees were used "
                             "to compute any reliable oob estimates.")
                    decision = predictions[k] \
                               / predictions[k].sum(axis=1)[:, np.newaxis]
                    self.oob_decision_function_.append(decision)

                    self.oob_score_ += np.mean(y[:, k] \
                                       == np.argmax(predictions[k], axis=1))

                if self.n_outputs_ == 1:
                    self.oob_decision_function_ = \
                        self.oob_decision_function_[0]

                self.oob_score_ /= self.n_outputs_

            else:
                # Regression:
                predictions = np.zeros((n_samples, self.n_outputs_))
                n_predictions = np.zeros((n_samples, self.n_outputs_))

                for estimator in self.estimators_:
                    mask = np.ones(n_samples, dtype=np.bool)
                    mask[estimator.indices_] = False

                    p_estimator = estimator.predict(X[mask, :])
                    if self.n_outputs_ == 1:
                        p_estimator = p_estimator[:, np.newaxis]

                    predictions[mask, :] += p_estimator
                    n_predictions[mask, :] += 1
                if (n_predictions == 0).any():
                    warn("Some inputs do not have OOB scores. "
                         "This probably means too few trees were used "
                         "to compute any reliable oob estimates.")
                    n_predictions[n_predictions == 0] = 1
                predictions /= n_predictions

                self.oob_prediction_ = predictions
                if self.n_outputs_ == 1:
                    self.oob_prediction_ = \
                        self.oob_prediction_.reshape((n_samples, ))

                self.oob_score_ = 0.0
                for k in xrange(self.n_outputs_):
                    self.oob_score_ += r2_score(y[:, k], predictions[:, k])
                self.oob_score_ /= self.n_outputs_

        # Sum the importances
        if self.compute_importances:
            self.feature_importances_ = \
                sum(tree.feature_importances_ for tree in self.estimators_) \
                / self.n_estimators

        return self


class ForestClassifier(BaseForest, ClassifierMixin):
    """Base class for forest of trees-based classifiers.

    Warning: This class should not be used directly. Use derived classes
    instead.
    """
    __metaclass__ = ABCMeta

    @abstractmethod
    def __init__(self, base_estimator,
                       n_estimators=10,
                       estimator_params=[],
                       bootstrap=False,
                       compute_importances=False,
                       oob_score=False,
                       n_jobs=1,
                       random_state=None,
                       verbose=0):
        super(ForestClassifier, self).__init__(
            base_estimator,
            n_estimators=n_estimators,
            estimator_params=estimator_params,
            bootstrap=bootstrap,
            compute_importances=compute_importances,
            oob_score=oob_score,
            n_jobs=n_jobs,
            random_state=random_state,
            verbose=verbose)

    def predict(self, X):
        """Predict class for X.

        The predicted class of an input sample is computed as the majority
        prediction of the trees in the forest.

        Parameters
        ----------
        X : array-like of shape = [n_samples, n_features]
            The input samples.

        Returns
        -------
        y : array of shape = [n_samples] or [n_samples, n_outputs]
            The predicted classes.
        """
        n_samples = len(X)

        P = self.predict_proba(X)
        if self.n_outputs_ == 1:
            P = [P]

        predictions = np.zeros((n_samples, self.n_outputs_))

        for k in xrange(self.n_outputs_):
            predictions[:, k] = self.classes_[k].take(np.argmax(P[k], axis=1),
                                                      axis=0)

        if self.n_outputs_ == 1:
            predictions = predictions.reshape((n_samples, ))

        return predictions

    def predict_proba(self, X):
        """Predict class probabilities for X.

        The predicted class probabilities of an input sample is computed as
        the mean predicted class probabilities of the trees in the forest.

        Parameters
        ----------
        X : array-like of shape = [n_samples, n_features]
            The input samples.

        Returns
        -------
        p : array of shape = [n_samples, n_classes], or a list of n_outputs
            such arrays if n_outputs > 1.
            The class probabilities of the input samples. Classes are
            ordered by arithmetical order.
        """
        # Check data
        if getattr(X, "dtype", None) != DTYPE or X.ndim != 2:
            X = array2d(X, dtype=DTYPE)

        # Assign chunk of trees to jobs
        n_jobs, n_trees, starts = _partition_trees(self)

        # Parallel loop
        all_p = Parallel(n_jobs=n_jobs)(
            delayed(_parallel_predict_proba)(
                self.estimators_[starts[i]:starts[i + 1]],
                X,
                self.n_classes_,
                self.n_outputs_)
            for i in xrange(n_jobs))

        # Reduce
        p = all_p[0]

        for j in xrange(1, len(all_p)):
            for k in xrange(self.n_outputs_):
                p[k] += all_p[j][k]

        for k in xrange(self.n_outputs_):
            p[k] /= self.n_estimators

        if self.n_outputs_ == 1:
            return p[0]

        else:
            return p

    def predict_log_proba(self, X):
        """Predict class log-probabilities for X.

        The predicted class log-probabilities of an input sample is computed as
        the mean predicted class log-probabilities of the trees in the forest.

        Parameters
        ----------
        X : array-like of shape = [n_samples, n_features]
            The input samples.

        Returns
        -------
        p : array of shape = [n_samples, n_classes], or a list of n_outputs
            such arrays if n_outputs > 1.
            The class log-probabilities of the input samples. Classes are
            ordered by arithmetical order.
        """
        proba = self.predict_proba(X)

        if self.n_outputs_ == 1:
            return np.log(proba)

        else:
            for k in xrange(self.n_outputs_):
                proba[k] = np.log(proba[k])

            return proba


class ForestRegressor(BaseForest, RegressorMixin):
    """Base class for forest of trees-based regressors.

    Warning: This class should not be used directly. Use derived classes
    instead.
    """
    __metaclass__ = ABCMeta

    @abstractmethod
    def __init__(self, base_estimator,
                       n_estimators=10,
                       estimator_params=[],
                       bootstrap=False,
                       compute_importances=False,
                       oob_score=False,
                       n_jobs=1,
                       random_state=None,
                       verbose=0):
        super(ForestRegressor, self).__init__(
            base_estimator,
            n_estimators=n_estimators,
            estimator_params=estimator_params,
            bootstrap=bootstrap,
            compute_importances=compute_importances,
            oob_score=oob_score,
            n_jobs=n_jobs,
            random_state=random_state,
            verbose=verbose)

    def predict(self, X):
        """Predict regression target for X.

        The predicted regression target of an input sample is computed as the
        mean predicted regression targets of the trees in the forest.

        Parameters
        ----------
        X : array-like of shape = [n_samples, n_features]
            The input samples.

        Returns
        -------
        y: array of shape = [n_samples] or [n_samples, n_outputs]
            The predicted values.
        """
        # Check data
        if getattr(X, "dtype", None) != DTYPE or X.ndim != 2:
            X = array2d(X, dtype=DTYPE)

        # Assign chunk of trees to jobs
        n_jobs, n_trees, starts = _partition_trees(self)

        # Parallel loop
        all_y_hat = Parallel(n_jobs=n_jobs)(
            delayed(_parallel_predict_regression)(
                self.estimators_[starts[i]:starts[i + 1]], X)
            for i in xrange(n_jobs))

        # Reduce
        y_hat = sum(all_y_hat) / self.n_estimators

        return y_hat


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
        all leaves are pure or until all leaves contain less than
        min_samples_split samples.
        Note: this parameter is tree-specific.

    min_samples_split : integer, optional (default=1)
        The minimum number of samples required to split an internal node.
        Note: this parameter is tree-specific.

    min_samples_leaf : integer, optional (default=1)
        The minimum number of samples in newly created leaves.  A split is
        discarded if after the split, one of the leaves would contain less then
        ``min_samples_leaf`` samples.
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

    compute_importances : boolean, optional (default=True)
        Whether feature importances are computed and stored into the
        ``feature_importances_`` attribute when calling fit.

    oob_score : bool
        Whether to use out-of-bag samples to estimate
        the generalization error.

    n_jobs : integer, optional (default=1)
        The number of jobs to run in parallel. If -1, then the number of jobs
        is set to the number of cores.

    random_state : int, RandomState instance or None, optional (default=None)
        If int, random_state is the seed used by the random number generator;
        If RandomState instance, random_state is the random number generator;
        If None, the random number generator is the RandomState instance used
        by `np.random`.

    verbose : int, optional (default=0)
        Controlls the verbosity of the tree building process.

    Attributes
    ----------
    `estimators_`: list of DecisionTreeClassifier
        The collection of fitted sub-estimators.

    `feature_importances_` : array, shape = [n_features]
        The feature importances (the higher, the more important the feature).

    `oob_score_` : float
        Score of the training dataset obtained using an out-of-bag estimate.

    `oob_decision_function_` : array, shape = [n_samples, n_classes]
        Decision function computed with out-of-bag estimate on the training
        set.


    References
    ----------

    .. [1] L. Breiman, "Random Forests", Machine Learning, 45(1), 5-32, 2001.

    See also
    --------
    DecisionTreeClassifier, ExtraTreesClassifier
    """
    def __init__(self, n_estimators=10,
                       criterion="gini",
                       max_depth=None,
                       min_samples_split=1,
                       min_samples_leaf=1,
                       min_density=0.1,
                       max_features="auto",
                       bootstrap=True,
                       compute_importances=False,
                       oob_score=False,
                       n_jobs=1,
                       random_state=None,
                       verbose=0):
        super(RandomForestClassifier, self).__init__(
            base_estimator=DecisionTreeClassifier(),
            n_estimators=n_estimators,
            estimator_params=("criterion", "max_depth", "min_samples_split",
                "min_samples_leaf", "min_density", "max_features",
                "random_state"),
            bootstrap=bootstrap,
            compute_importances=compute_importances,
            oob_score=oob_score,
            n_jobs=n_jobs,
            random_state=random_state,
            verbose=verbose)

        self.criterion = criterion
        self.max_depth = max_depth
        self.min_samples_split = min_samples_split
        self.min_samples_leaf = min_samples_leaf
        self.min_density = min_density
        self.max_features = max_features


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
        all leaves are pure or until all leaves contain less than
        min_samples_split samples.
        Note: this parameter is tree-specific.

    min_samples_split : integer, optional (default=1)
        The minimum number of samples required to split an internal node.
        Note: this parameter is tree-specific.

    min_samples_leaf : integer, optional (default=1)
        The minimum number of samples in newly created leaves.  A split is
        discarded if after the split, one of the leaves would contain less then
        ``min_samples_leaf`` samples.
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

    compute_importances : boolean, optional (default=True)
        Whether feature importances are computed and stored into the
        ``feature_importances_`` attribute when calling fit.

    oob_score : bool
        whether to use out-of-bag samples to estimate
        the generalization error.

    n_jobs : integer, optional (default=1)
        The number of jobs to run in parallel. If -1, then the number of jobs
        is set to the number of cores.

    random_state : int, RandomState instance or None, optional (default=None)
        If int, random_state is the seed used by the random number generator;
        If RandomState instance, random_state is the random number generator;
        If None, the random number generator is the RandomState instance used
        by `np.random`.

    verbose : int, optional (default=0)
        Controlls the verbosity of the tree building process.

    Attributes
    ----------
    `estimators_`: list of DecisionTreeRegressor
        The collection of fitted sub-estimators.

    `feature_importances_` : array of shape = [n_features]
        The feature mportances (the higher, the more important the feature).

    `oob_score_` : float
        Score of the training dataset obtained using an out-of-bag estimate.

    `oob_prediction_` : array, shape = [n_samples]
        Prediction computed with out-of-bag estimate on the training set.



    References
    ----------

    .. [1] L. Breiman, "Random Forests", Machine Learning, 45(1), 5-32, 2001.

    See also
    --------
    DecisionTreeRegressor, ExtraTreesRegressor
    """
    def __init__(self, n_estimators=10,
                       criterion="mse",
                       max_depth=None,
                       min_samples_split=1,
                       min_samples_leaf=1,
                       min_density=0.1,
                       max_features="auto",
                       bootstrap=True,
                       compute_importances=False,
                       oob_score=False,
                       n_jobs=1,
                       random_state=None,
                       verbose=0):
        super(RandomForestRegressor, self).__init__(
            base_estimator=DecisionTreeRegressor(),
            n_estimators=n_estimators,
            estimator_params=("criterion", "max_depth", "min_samples_split",
                "min_samples_leaf", "min_density", "max_features",
                "random_state"),
            bootstrap=bootstrap,
            compute_importances=compute_importances,
            oob_score=oob_score,
            n_jobs=n_jobs,
            random_state=random_state,
            verbose=verbose)

        self.criterion = criterion
        self.max_depth = max_depth
        self.min_samples_split = min_samples_split
        self.min_samples_leaf = min_samples_leaf
        self.min_density = min_density
        self.max_features = max_features


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
        all leaves are pure or until all leaves contain less than
        min_samples_split samples.
        Note: this parameter is tree-specific.

    min_samples_split : integer, optional (default=1)
        The minimum number of samples required to split an internal node.
        Note: this parameter is tree-specific.

    min_samples_leaf : integer, optional (default=1)
        The minimum number of samples in newly created leaves.  A split is
        discarded if after the split, one of the leaves would contain less then
        ``min_samples_leaf`` samples.
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

    compute_importances : boolean, optional (default=True)
        Whether feature importances are computed and stored into the
        ``feature_importances_`` attribute when calling fit.

    oob_score : bool
        Whether to use out-of-bag samples to estimate
        the generalization error.

    n_jobs : integer, optional (default=1)
        The number of jobs to run in parallel. If -1, then the number of jobs
        is set to the number of cores.

    random_state : int, RandomState instance or None, optional (default=None)
        If int, random_state is the seed used by the random number generator;
        If RandomState instance, random_state is the random number generator;
        If None, the random number generator is the RandomState instance used
        by `np.random`.

    verbose : int, optional (default=0)
        Controlls the verbosity of the tree building process.

    Attributes
    ----------
    `estimators_`: list of DecisionTreeClassifier
        The collection of fitted sub-estimators.

    `feature_importances_` : array of shape = [n_features]
        The feature mportances (the higher, the more important the feature).

    `oob_score_` : float
        Score of the training dataset obtained using an out-of-bag estimate.

    `oob_decision_function_` : array, shape = [n_samples, n_classes]
        Decision function computed with out-of-bag estimate on the training
        set.

    References
    ----------

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
                       min_samples_split=1,
                       min_samples_leaf=1,
                       min_density=0.1,
                       max_features="auto",
                       bootstrap=False,
                       compute_importances=False,
                       oob_score=False,
                       n_jobs=1,
                       random_state=None,
                       verbose=0):
        super(ExtraTreesClassifier, self).__init__(
            base_estimator=ExtraTreeClassifier(),
            n_estimators=n_estimators,
            estimator_params=("criterion", "max_depth", "min_samples_split",
                "min_samples_leaf", "min_density", "max_features",
                "random_state"),
            bootstrap=bootstrap,
            compute_importances=compute_importances,
            oob_score=oob_score,
            n_jobs=n_jobs,
            random_state=random_state,
            verbose=verbose)

        self.criterion = criterion
        self.max_depth = max_depth
        self.min_samples_split = min_samples_split
        self.min_samples_leaf = min_samples_leaf
        self.min_density = min_density
        self.max_features = max_features


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
        all leaves are pure or until all leaves contain less than
        min_samples_split samples.
        Note: this parameter is tree-specific.

    min_samples_split : integer, optional (default=1)
        The minimum number of samples required to split an internal node.
        Note: this parameter is tree-specific.

    min_samples_leaf : integer, optional (default=1)
        The minimum number of samples in newly created leaves.  A split is
        discarded if after the split, one of the leaves would contain less then
        ``min_samples_leaf`` samples.
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

    compute_importances : boolean, optional (default=True)
        Whether feature importances are computed and stored into the
        ``feature_importances_`` attribute when calling fit.

    oob_score : bool
        Whether to use out-of-bag samples to estimate
        the generalization error.

    n_jobs : integer, optional (default=1)
        The number of jobs to run in parallel. If -1, then the number of jobs
        is set to the number of cores.

    random_state : int, RandomState instance or None, optional (default=None)
        If int, random_state is the seed used by the random number generator;
        If RandomState instance, random_state is the random number generator;
        If None, the random number generator is the RandomState instance used
        by `np.random`.

    verbose : int, optional (default=0)
        Controlls the verbosity of the tree building process.

    Attributes
    ----------
    `estimators_`: list of DecisionTreeRegressor
        The collection of fitted sub-estimators.

    `feature_importances_` : array of shape = [n_features]
        The feature mportances (the higher, the more important the feature).

    `oob_score_` : float
        Score of the training dataset obtained using an out-of-bag estimate.

    `oob_prediction_` : array, shape = [n_samples]
        Prediction computed with out-of-bag estimate on the training set.

    References
    ----------

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
                       min_samples_split=1,
                       min_samples_leaf=1,
                       min_density=0.1,
                       max_features="auto",
                       bootstrap=False,
                       compute_importances=False,
                       oob_score=False,
                       n_jobs=1,
                       random_state=None,
                       verbose=0):
        super(ExtraTreesRegressor, self).__init__(
            base_estimator=ExtraTreeRegressor(),
            n_estimators=n_estimators,
            estimator_params=("criterion", "max_depth", "min_samples_split",
                "min_samples_leaf", "min_density", "max_features",
                "random_state"),
            bootstrap=bootstrap,
            compute_importances=compute_importances,
            oob_score=oob_score,
            n_jobs=n_jobs,
            random_state=random_state,
            verbose=verbose)

        self.criterion = criterion
        self.max_depth = max_depth
        self.min_samples_split = min_samples_split
        self.min_samples_leaf = min_samples_leaf
        self.min_density = min_density
        self.max_features = max_features
