"""Forest of trees-based ensemble methods

Those methods include random forests and extremely randomized trees.

The module structure is the following:

- The ``Forest`` base class implements a common ``fit`` method for all
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

import numpy as np

from ..base import ClassifierMixin, RegressorMixin
from ..externals.joblib import Parallel, delayed
from ..tree import DecisionTreeClassifier, DecisionTreeRegressor, \
                   ExtraTreeClassifier, ExtraTreeRegressor
from ..utils import check_random_state

from .base import BaseEnsemble

__all__ = ["RandomForestClassifier",
           "RandomForestRegressor",
           "ExtraTreesClassifier",
           "ExtraTreesRegressor"]


def _parallel_build_tree(forest, X, y, sample_mask, X_argsorted):
    tree = forest._make_estimator(append=False)

    if forest.bootstrap:
        n_samples = X.shape[0]
        indices = forest.random_state.randint(0, n_samples, n_samples)
        tree.fit(X[indices], y[indices],
                 sample_mask=sample_mask, X_argsorted=X_argsorted)

    else:
        tree.fit(X, y,
                 sample_mask=sample_mask, X_argsorted=X_argsorted)

    return tree

def _parallel_predict_proba(tree, X, n_classes):
    p = np.zeros((X.shape[0], n_classes))

    if n_classes == tree.n_classes_:
        p += tree.predict_proba(X)

    else:
        proba = tree.predict_proba(X)

        for j, c in enumerate(tree.classes_):
            p[:, c] += proba[:, j]

    return p

def _parallel_predict_regr(tree, X):
    return tree.predict(X)

def _parallel_compute_importances(tree):
    return tree.feature_importances()


class Forest(BaseEnsemble):
    """Base class for forests of trees.

    Warning: This class should not be used directly. Use derived classes
    instead.
    """
    def __init__(self, base_estimator,
                       n_estimators=10,
                       estimator_params=[],
                       bootstrap=False,
                       n_jobs=1,
                       random_state=None):
        super(Forest, self).__init__(
            base_estimator=base_estimator,
            n_estimators=n_estimators,
            estimator_params=estimator_params)

        self.bootstrap = bootstrap
        self.n_jobs = n_jobs
        self.random_state = check_random_state(random_state)

    def fit(self, X, y):
        """Build a forest of trees from the training set (X, y).

        Parameters
        ----------
        X : array-like of shape = [n_samples, n_features]
            The training input samples.

        y : array-like, shape = [n_samples]
            The target values (integers that correspond to classes in
            classification, real numbers in regression).

        Returns
        -------
        self : object
            Returns self.
        """
        # Build the forest
        X = np.atleast_2d(X)
        y = np.atleast_1d(y)

        if self.bootstrap:
            sample_mask = None
            X_argsorted = None

        else:
            sample_mask = np.ones((X.shape[0],), dtype=np.bool)
            X_argsorted = np.asfortranarray(
                np.argsort(X.T, axis=1).astype(np.int32).T)

        if isinstance(self.base_estimator, ClassifierMixin):
            self.classes_ = np.unique(y)
            self.n_classes_ = len(self.classes_)
            y = np.searchsorted(self.classes_, y)

        self.estimators_ = Parallel(n_jobs=self.n_jobs)(
            delayed(_parallel_build_tree)(self, X, y, sample_mask, X_argsorted)
                for i in xrange(self.n_estimators))

        return self

    def feature_importances(self):
        """Compute the mean feature importances over the trees in the forest.

        Returns
        -------
        importances : array of shape = [n_features]
            The feature importances.
        """
        all_importances = Parallel(n_jobs=self.n_jobs)(
            delayed(_parallel_compute_importances)(self[i])
                for i in xrange(self.n_estimators))

        importances = sum(all_importances) / self.n_estimators

        return importances


class ForestClassifier(Forest, ClassifierMixin):
    """Base class for forest of trees-based classifiers.

    Warning: This class should not be used directly. Use derived classes
    instead.
    """
    def __init__(self, base_estimator,
                       n_estimators=10,
                       estimator_params=[],
                       bootstrap=False,
                       n_jobs=1,
                       random_state=None):
        super(ForestClassifier, self).__init__(
            base_estimator,
            n_estimators=n_estimators,
            estimator_params=estimator_params,
            bootstrap=bootstrap,
            n_jobs=n_jobs,
            random_state=random_state)

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
        y : array of shape = [n_samples]
            The predicted classes.
        """
        return self.classes_.take(
            np.argmax(self.predict_proba(X), axis=1),  axis=0)

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
        p : array of shape = [n_samples]
            The class probabilities of the input samples. Classes are
            ordered by arithmetical order.
        """
        X = np.atleast_2d(X)

        all_p = Parallel(n_jobs=self.n_jobs)(
            delayed(_parallel_predict_proba)(self[i], X, self.n_classes_)
                for i in xrange(self.n_estimators))

        p = sum(all_p) / self.n_estimators

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
        p : array of shape = [n_samples]
            The class log-probabilities of the input samples. Classes are
            ordered by arithmetical order.
        """
        return np.log(self.predict_proba(X))


class ForestRegressor(Forest, RegressorMixin):
    """Base class for forest of trees-based regressors.

    Warning: This class should not be used directly. Use derived classes
    instead.
    """
    def __init__(self, base_estimator,
                       n_estimators=10,
                       estimator_params=[],
                       bootstrap=False,
                       n_jobs=1,
                       random_state=None):
        super(ForestRegressor, self).__init__(
            base_estimator,
            n_estimators=n_estimators,
            estimator_params=estimator_params,
            bootstrap=bootstrap,
            n_jobs=n_jobs,
            random_state=random_state)

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
        y: array of shape = [n_samples]
            The predicted values.
        """
        X = np.atleast_2d(X)

        all_y_hat = Parallel(n_jobs=self.n_jobs)(
            delayed(_parallel_predict_regr)(self[i], X)
                for i in xrange(self.n_estimators))

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

    max_depth : integer or None, optional (default=10)
        The maximum depth of the tree. If None, then nodes are expanded until
        all leaves are pure or until all leaves contain less than min_split
        samples.

    min_split : integer, optional (default=1)
        The minimum number of samples required to split an internal node.

    min_density : float, optional (default=0.1)
        The minimum density of the `sample_mask` (i.e. the fraction of samples
        in the mask). If the density falls below this threshold the mask is
        recomputed and the input data is packed which results in data copying.
        If `min_density` equals to one, the partitions are always represented
        as copies of the original data. Otherwise, partitions are represented
        as bit masks (a.k.a. sample masks).

    max_features : int or None, optional (default=None)
        The number of features to consider when looking for the best split.
        If None, all features are considered, otherwise max_features are chosen
        at random.

    bootstrap : boolean, optional (default=True)
        Whether bootstrap samples are used when building trees.

    random_state : int, RandomState instance or None, optional (default=None)
        If int, random_state is the seed used by the random number generator;
        If RandomState instance, random_state is the random number generator;
        If None, the random number generator is the RandomState instance used
        by `np.random`.

    See also
    --------
    RandomForestRegressor, ExtraTreesClassifier, ExtraTreesRegressor

    References
    ----------
    .. [1] L. Breiman, "Random Forests", Machine Learning, 45(1), 5-32, 2001.
    """
    def __init__(self, n_estimators=10,
                       criterion="gini",
                       max_depth=10,
                       min_split=1,
                       min_density=0.1,
                       max_features=None,
                       bootstrap=True,
                       n_jobs=1,
                       random_state=None):
        super(RandomForestClassifier, self).__init__(
            base_estimator=DecisionTreeClassifier(),
            n_estimators=n_estimators,
            estimator_params=("criterion", "max_depth", "min_split",
                              "min_density", "max_features", "random_state"),
            bootstrap=bootstrap,
            n_jobs=n_jobs,
            random_state=random_state)

        self.criterion = criterion
        self.max_depth = max_depth
        self.min_split = min_split
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

    max_depth : integer or None, optional (default=10)
        The maximum depth of the tree. If None, then nodes are expanded until
        all leaves are pure or until all leaves contain less than min_split
        samples.

    min_split : integer, optional (default=1)
        The minimum number of samples required to split an internal node.

    min_density : float, optional (default=0.1)
        The minimum density of the `sample_mask` (i.e. the fraction of samples
        in the mask). If the density falls below this threshold the mask is
        recomputed and the input data is packed which results in data copying.
        If `min_density` equals to one, the partitions are always represented
        as copies of the original data. Otherwise, partitions are represented
        as bit masks (a.k.a. sample masks).

    max_features : int or None, optional (default=None)
        The number of features to consider when looking for the best split.
        If None, all features are considered, otherwise max_features are chosen
        at random.

    bootstrap : boolean, optional (default=True)
        Whether bootstrap samples are used when building trees.

    random_state : int, RandomState instance or None, optional (default=None)
        If int, random_state is the seed used by the random number generator;
        If RandomState instance, random_state is the random number generator;
        If None, the random number generator is the RandomState instance used
        by `np.random`.

    See also
    --------
    RandomForestClassifier, ExtraTreesClassifier, ExtraTreesRegressor

    References
    ----------
    .. [1] L. Breiman, "Random Forests", Machine Learning, 45(1), 5-32, 2001.
    """
    def __init__(self, n_estimators=10,
                       criterion="mse",
                       max_depth=10,
                       min_split=1,
                       min_density=0.1,
                       max_features=None,
                       bootstrap=True,
                       n_jobs=1,
                       random_state=None):
        super(RandomForestRegressor, self).__init__(
            base_estimator=DecisionTreeRegressor(),
            n_estimators=n_estimators,
            estimator_params=("criterion", "max_depth", "min_split",
                              "min_density", "max_features", "random_state"),
            bootstrap=bootstrap,
            n_jobs=n_jobs,
            random_state=random_state)

        self.criterion = criterion
        self.max_depth = max_depth
        self.min_split = min_split
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

    max_depth : integer or None, optional (default=10)
        The maximum depth of the tree. If None, then nodes are expanded until
        all leaves are pure or until all leaves contain less than min_split
        samples.

    min_split : integer, optional (default=1)
        The minimum number of samples required to split an internal node.

    min_density : float, optional (default=0.1)
        The minimum density of the `sample_mask` (i.e. the fraction of samples
        in the mask). If the density falls below this threshold the mask is
        recomputed and the input data is packed which results in data copying.
        If `min_density` equals to one, the partitions are always represented
        as copies of the original data. Otherwise, partitions are represented
        as bit masks (a.k.a. sample masks).

    max_features : int or None, optional (default=None)
        The number of features to consider when looking for the best split.
        If None, all features are considered, otherwise max_features are chosen
        at random.

    bootstrap : boolean, optional (default=False)
        Whether bootstrap samples are used when building trees.

    random_state : int, RandomState instance or None, optional (default=None)
        If int, random_state is the seed used by the random number generator;
        If RandomState instance, random_state is the random number generator;
        If None, the random number generator is the RandomState instance used
        by `np.random`.

    See also
    --------
    ExtraTreesRegressor, RandomForestClassifier, RandomForestRegressor

    References
    ----------
    .. [1] P. Geurts, D. Ernst., and L. Wehenkel, "Extremely randomized trees",
           Machine Learning, 63(1), 3-42, 2006.
    """
    def __init__(self, n_estimators=10,
                       criterion="gini",
                       max_depth=10,
                       min_split=1,
                       min_density=0.1,
                       max_features=None,
                       bootstrap=False,
                       n_jobs=1,
                       random_state=None):
        super(ExtraTreesClassifier, self).__init__(
            base_estimator=ExtraTreeClassifier(),
            n_estimators=n_estimators,
            estimator_params=("criterion", "max_depth", "min_split",
                              "min_density", "max_features", "random_state"),
            bootstrap=bootstrap,
            n_jobs=n_jobs,
            random_state=random_state)

        self.criterion = criterion
        self.max_depth = max_depth
        self.min_split = min_split
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

    max_depth : integer or None, optional (default=10)
        The maximum depth of the tree. If None, then nodes are expanded until
        all leaves are pure or until all leaves contain less than min_split
        samples.

    min_split : integer, optional (default=1)
        The minimum number of samples required to split an internal node.

    min_density : float, optional (default=0.1)
        The minimum density of the `sample_mask` (i.e. the fraction of samples
        in the mask). If the density falls below this threshold the mask is
        recomputed and the input data is packed which results in data copying.
        If `min_density` equals to one, the partitions are always represented
        as copies of the original data. Otherwise, partitions are represented
        as bit masks (a.k.a. sample masks).

    max_features : int or None, optional (default=None)
        The number of features to consider when looking for the best split.
        If None, all features are considered, otherwise max_features are chosen
        at random.

    bootstrap : boolean, optional (default=False)
        Whether bootstrap samples are used when building trees.

    random_state : int, RandomState instance or None, optional (default=None)
        If int, random_state is the seed used by the random number generator;
        If RandomState instance, random_state is the random number generator;
        If None, the random number generator is the RandomState instance used
        by `np.random`.

    See also
    --------
    ExtraTreesRegressor, RandomForestClassifier, RandomForestRegressor

    References
    ----------
    .. [1] P. Geurts, D. Ernst., and L. Wehenkel, "Extremely randomized trees",
           Machine Learning, 63(1), 3-42, 2006.
    """
    def __init__(self, n_estimators=10,
                       criterion="mse",
                       max_depth=10,
                       min_split=1,
                       min_density=0.1,
                       max_features=None,
                       bootstrap=False,
                       n_jobs=1,
                       random_state=None):
        super(ExtraTreesRegressor, self).__init__(
            base_estimator=ExtraTreeRegressor(),
            n_estimators=n_estimators,
            estimator_params=("criterion", "max_depth", "min_split",
                              "min_density", "max_features", "random_state"),
            bootstrap=bootstrap,
            n_jobs=n_jobs,
            random_state=random_state)

        self.criterion = criterion
        self.max_depth = max_depth
        self.min_split = min_split
        self.min_density = min_density
        self.max_features = max_features
