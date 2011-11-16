"""
This module gathers forest of trees-based ensemble methods, including random
forests and extra-trees.
"""

# Authors: Gilles Louppe, Brian Holt
# License: BSD 3

import numpy as np

from ..base import BaseEstimator, ClassifierMixin, RegressorMixin
from ..tree import DecisionTreeClassifier, DecisionTreeRegressor, \
                   ExtraTreeClassifier, ExtraTreeRegressor
from ..utils import check_random_state

from .base import BaseEnsemble

__all__ = ["RandomForestClassifier",
           "RandomForestRegressor",
           "ExtraTreesClassifier",
           "ExtraTreesRegressor"]


class Forest(BaseEnsemble):
    """Base class for forests of trees.

    Warning: This class should not be used directly. Use derived classes
    instead.
    """
    def __init__(self, base_tree,
                       n_trees=10,
                       bootstrap=False,
                       random_state=None,
                       **tree_params):
        super(Forest, self).__init__(
            base_estimator=base_tree,
            n_estimators=n_trees,
            **tree_params)

        self.bootstrap = bootstrap
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

        Return
        ------
        self : object
            Returns self.
        """
        # Build the forest
        X = np.atleast_2d(X)
        y = np.atleast_1d(y)

        if isinstance(self.base_estimator, ClassifierMixin):
            self.classes = np.unique(y)
            self.n_classes = len(self.classes)
            y = np.searchsorted(self.classes, y)

        for tree in self.estimators:
            tree.set_params(random_state=self.random_state)

            if self.bootstrap:
                n_samples = X.shape[0]
                indices = self.random_state.randint(0, n_samples, n_samples)

                X = X[indices]
                y = y[indices]

            tree.fit(X, y)

        return self


class ForestClassifier(Forest, ClassifierMixin):
    """Base class for forest of trees-based classifiers.

    Warning: This class should not be used directly. Use derived classes
    instead.
    """
    def __init__(self, base_tree,
                       n_trees=10,
                       bootstrap=False,
                       random_state=None,
                       **tree_params):
        super(ForestClassifier, self).__init__(
            base_tree,
            n_trees,
            bootstrap=bootstrap,
            random_state=random_state,
            **tree_params)

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
        predictions : array of shape = [n_samples]
            The predicted classes.
        """
        return self.classes.take(
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
        p = np.zeros((X.shape[0], self.n_classes))

        for tree in self.estimators:
            p += tree.predict_proba(X)

        p /= self.n_estimators

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
    def __init__(self, base_tree,
                       n_trees=10,
                       bootstrap=False,
                       random_state=None,
                       **tree_params):
        super(ForestRegressor, self).__init__(
            base_tree,
            n_trees,
            bootstrap=bootstrap,
            random_state=random_state,
            **tree_params)

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
        predictions : array of shape = [n_samples]
            The predicted values.
        """
        X = np.atleast_2d(X)
        y_hat = np.zeros(X.shape[0])

        for tree in self.estimators:
            y_hat += tree.predict(X)

        y_hat /= self.n_estimators

        return y_hat


class RandomForestClassifier(ForestClassifier):
    """A random forest classifier [1].

    Parameters
    ----------
    base_tree : object, optional (default=None)
        The base tree from which the forest is built. If None, a
        `DecisionTreeClassifier` with parameters defined from **tree_params
        is used.

    n_trees : integer, optional (default=10)
        The number of trees in the forest.

    bootstrap : boolean, optional (default=True)
        Whether bootstrap samples are used when building trees.

    random_state : int, RandomState instance or None, optional (default=None)
        If int, random_state is the seed used by the random number generator;
        If RandomState instance, random_state is the random number generator;
        If None, the random number generator is the RandomState instance used
        by `np.random`.

    **tree_params : key-words parameters
        The parameters to set when instantiating the underlying base tree. If
        none are given, default parameters are used.

    Attributes
    ----------
    base_tree : object
        The underlying tree that is used to generate the forest.

    Notes
    -----
    When using grid search to optimize the parameters, use the nested object
    syntax to set the parameters of the underlying trees (e.g.,
    `base_tree__max_depth`).

    See also
    --------
    RandomForestRegressor, ExtraTreesClassifier, ExtraTreesRegressor

    References
    ----------
    .. [1] L. Breiman, "Random Forests", Machine Learning, 45(1), 5-32, 2001.
    """
    def __init__(self, base_tree=None,
                       n_trees=10,
                       bootstrap=True,
                       random_state=None,
                       **tree_params):
        super(RandomForestClassifier, self).__init__(
            base_tree if base_tree is not None \
                      else DecisionTreeClassifier(),
            n_trees,
            bootstrap=bootstrap,
            random_state=random_state,
            **tree_params)


class RandomForestRegressor(ForestRegressor):
    """A random forest regressor [1].

    Parameters
    ----------
    base_tree : object, optional (default=None)
        The base tree from which the forest is built. If None, a
        `DecisionTreeRegressor` with parameters defined from **tree_params
        is used.

    n_trees : integer, optional (default=10)
        The number of trees in the forest.

    bootstrap : boolean, optional (default=True)
        Whether bootstrap samples are used when building trees.

    random_state : int, RandomState instance or None, optional (default=None)
        If int, random_state is the seed used by the random number generator;
        If RandomState instance, random_state is the random number generator;
        If None, the random number generator is the RandomState instance used
        by `np.random`.

    **tree_params : key-words parameters
        The parameters to set when instantiating the underlying base tree. If
        none are given, default parameters are used.

    Attributes
    ----------
    base_tree : object
        The underlying tree that is used to generate the forest.

    Notes
    -----
    When using grid search to optimize the parameters, use the nested object
    syntax to set the parameters of the underlying trees (e.g.,
    `base_tree__max_depth`).

    See also
    --------
    RandomForestClassifier, ExtraTreesClassifier, ExtraTreesRegressor

    References
    ----------
    .. [1] L. Breiman, "Random Forests", Machine Learning, 45(1), 5-32, 2001.
    """
    def __init__(self, base_tree=None,
                       n_trees=10,
                       bootstrap=True,
                       random_state=None,
                       **tree_params):
        super(RandomForestRegressor, self).__init__(
            base_tree if base_tree is not None \
                      else DecisionTreeRegressor(),
            n_trees,
            bootstrap=bootstrap,
            random_state=random_state,
            **tree_params)


class ExtraTreesClassifier(ForestClassifier):
    """An extra-trees classifier [1].

    Parameters
    ----------
    base_tree : object, optional (default=None)
        The base tree from which the forest is built. If None, an
        `ExtraTreeClassifier` with parameters defined from **tree_params
        is used.

    n_trees : integer, optional (default=10)
        The number of trees in the forest.

    bootstrap : boolean, optional (default=True)
        Whether bootstrap samples are used when building trees.

    random_state : int, RandomState instance or None, optional (default=None)
        If int, random_state is the seed used by the random number generator;
        If RandomState instance, random_state is the random number generator;
        If None, the random number generator is the RandomState instance used
        by `np.random`.

    **tree_params : key-words parameters
        The parameters to set when instantiating the underlying base tree. If
        none are given, default parameters are used.

    Attributes
    ----------
    base_tree : object
        The underlying tree that is used to generate the forest.

    Notes
    -----
    When using grid search to optimize the parameters, use the nested object
    syntax to set the parameters of the underlying trees (e.g.,
    `base_tree__max_depth`).

    See also
    --------
    ExtraTreesRegressor, RandomForestClassifier, RandomForestRegressor

    References
    ----------
    .. [1] P. Geurts, D. Ernst., and L. Wehenkel, "Extremely randomized trees",
           Machine Learning, 63(1), 3-42, 2006.
    """
    def __init__(self, base_tree=None,
                       n_trees=10,
                       bootstrap=False,
                       random_state=None,
                       **tree_params):
        super(ExtraTreesClassifier, self).__init__(
            base_tree if base_tree is not None \
                      else ExtraTreeClassifier(),
            n_trees,
            bootstrap=bootstrap,
            random_state=random_state,
            **tree_params)


class ExtraTreesRegressor(ForestRegressor):
    """An extra-trees regressor [1].

    Parameters
    ----------
    base_tree : object, optional (default=None)
        The base tree from which the forest is built. If None, an
        `ExtraTreeRegressor` with parameters defined from **tree_params
        is used.

    n_trees : integer, optional (default=10)
        The number of trees in the forest.

    bootstrap : boolean, optional (default=True)
        Whether bootstrap samples are used when building trees.

    random_state : int, RandomState instance or None, optional (default=None)
        If int, random_state is the seed used by the random number generator;
        If RandomState instance, random_state is the random number generator;
        If None, the random number generator is the RandomState instance used
        by `np.random`.

    **tree_params : key-words parameters
        The parameters to set when instantiating the underlying base tree. If
        none are given, default parameters are used.

    Attributes
    ----------
    base_tree : object
        The underlying tree that is used to generate the forest.

    Notes
    -----
    When using grid search to optimize the parameters, use the nested object
    syntax to set the parameters of the underlying trees (e.g.,
    `base_tree__max_depth`).

    See also
    --------
    ExtraTreesRegressor, RandomForestClassifier, RandomForestRegressor

    References
    ----------
    .. [1] P. Geurts, D. Ernst., and L. Wehenkel, "Extremely randomized trees",
           Machine Learning, 63(1), 3-42, 2006.
    """
    def __init__(self, base_tree=None,
                       n_trees=10,
                       bootstrap=False,
                       random_state=None,
                       **tree_params):
        super(ExtraTreesRegressor, self).__init__(
            base_tree if base_tree is not None \
                      else ExtraTreeRegressor(),
            n_trees,
            bootstrap=bootstrap,
            random_state=random_state,
            **tree_params)
