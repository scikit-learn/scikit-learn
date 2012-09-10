"""
This module gathers tree-based methods, including decision, regression and
randomized trees. Single and multi-output problems are both handled.
"""

# Code is originally adapted from MILK: Machine Learning Toolkit
# Copyright (C) 2008-2011, Luis Pedro Coelho <luis@luispedro.org>
# License: MIT. See COPYING.MIT file in the milk distribution

# Authors: Brian Holt, Peter Prettenhofer, Satrajit Ghosh, Gilles Louppe
# License: BSD3

from __future__ import division
import numpy as np
from abc import ABCMeta, abstractmethod

from ..base import BaseEstimator, ClassifierMixin, RegressorMixin
from ..feature_selection.selector_mixin import SelectorMixin
from ..utils import array2d, check_random_state

from . import _tree


__all__ = ["DecisionTreeClassifier",
           "DecisionTreeRegressor",
           "ExtraTreeClassifier",
           "ExtraTreeRegressor"]

DTYPE = _tree.DTYPE
DOUBLE = _tree.DOUBLE

CLASSIFICATION = {
    "gini": _tree.Gini,
    "entropy": _tree.Entropy,
}

REGRESSION = {
    "mse": _tree.MSE,
}


def export_graphviz(decision_tree, out_file=None, feature_names=None):
    """Export a decision tree in DOT format.

    This function generates a GraphViz representation of the decision tree,
    which is then written into `out_file`. Once exported, graphical renderings
    can be generated using, for example::

        $ dot -Tps tree.dot -o tree.ps      (PostScript format)
        $ dot -Tpng tree.dot -o tree.png    (PNG format)

    Parameters
    ----------
    decision_tree : decision tree classifier
        The decision tree to be exported to graphviz.

    out : file object or string, optional (default=None)
        Handle or name of the output file.

    feature_names : list of strings, optional (default=None)
        Names of each of the features.

    Returns
    -------
    out_file : file object
        The file object to which the tree was exported.  The user is
        expected to `close()` this object when done with it.

    Examples
    --------
    >>> from sklearn.datasets import load_iris
    >>> from sklearn import tree

    >>> clf = tree.DecisionTreeClassifier()
    >>> iris = load_iris()

    >>> clf = clf.fit(iris.data, iris.target)
    >>> import tempfile
    >>> out_file = tree.export_graphviz(clf, out_file=tempfile.TemporaryFile())
    >>> out_file.close()
    """
    def node_to_str(tree, node_id):
        value = tree.value[node_id]
        if tree.n_outputs == 1:
            value = value[0, :]

        if tree.children_left[node_id] == _tree.TREE_LEAF:
            return "error = %.4f\\nsamples = %s\\nvalue = %s" \
                   % (tree.init_error[node_id],
                      tree.n_samples[node_id],
                      value)
        else:
            if feature_names is not None:
                feature = feature_names[tree.feature[node_id]]
            else:
                feature = "X[%s]" % tree.feature[node_id]

            return "%s <= %.4f\\nerror = %s\\nsamples = %s\\nvalue = %s" \
                   % (feature,
                      tree.threshold[node_id],
                      tree.init_error[node_id],
                      tree.n_samples[node_id],
                      value)

    def recurse(tree, node_id, parent=None):
        if node_id == _tree.TREE_LEAF:
            raise ValueError("Invalid node_id %s" % _tree.TREE_LEAF)

        left_child = tree.children_left[node_id]
        right_child = tree.children_right[node_id]

        # Add node with description
        out_file.write('%d [label="%s", shape="box"] ;\n' %
                (node_id, node_to_str(tree, node_id)))

        if parent is not None:
            # Add edge to parent
            out_file.write('%d -> %d ;\n' % (parent, node_id))

        if left_child != _tree.TREE_LEAF:  # and right_child != _tree.TREE_LEAF
            recurse(tree, left_child, node_id)
            recurse(tree, right_child, node_id)

    if out_file is None:
        out_file = open("tree.dot", "w")
    elif isinstance(out_file, basestring):
        out_file = open(out_file, "w")

    out_file.write("digraph Tree {\n")
    if isinstance(decision_tree, _tree.Tree):
        recurse(decision_tree, 0)
    else:
        recurse(decision_tree.tree_, 0)
    out_file.write("}")

    return out_file


class BaseDecisionTree(BaseEstimator, SelectorMixin):
    """Base class for decision trees.

    Warning: This class should not be used directly.
    Use derived classes instead.
    """
    __metaclass__ = ABCMeta

    @abstractmethod
    def __init__(self, criterion,
                       max_depth,
                       min_samples_split,
                       min_samples_leaf,
                       min_density,
                       max_features,
                       compute_importances,
                       random_state):
        self.criterion = criterion
        self.max_depth = max_depth
        self.min_samples_split = min_samples_split
        self.min_samples_leaf = min_samples_leaf
        self.min_density = min_density
        self.max_features = max_features
        self.compute_importances = compute_importances
        self.random_state = random_state

        self.n_features_ = None
        self.n_outputs_ = None
        self.classes_ = None
        self.n_classes_ = None
        self.find_split_ = _tree.TREE_SPLIT_BEST

        self.tree_ = None
        self.feature_importances_ = None

    def fit(self, X, y, sample_mask=None, X_argsorted=None):
        """Build a decision tree from the training set (X, y).

        Parameters
        ----------
        X : array-like of shape = [n_samples, n_features]
            The training input samples. Use ``dtype=np.float32``
            and ``order='F'`` for maximum efficiency.

        y : array-like, shape = [n_samples] or [n_samples, n_outputs]
            The target values (integers that correspond to classes in
            classification, real numbers in regression).
            Use ``dtype=np.float64`` and ``order='C'`` for maximum
            efficiency.

        Returns
        -------
        self : object
            Returns self.
        """
        self.random_state = check_random_state(self.random_state)

        # set min_samples_split sensibly
        self.min_samples_split = max(self.min_samples_split,
                                     2 * self.min_samples_leaf)

        # Convert data
        if getattr(X, "dtype", None) != DTYPE or \
           X.ndim != 2 or not X.flags.fortran:
            X = array2d(X, dtype=DTYPE, order="F")

        n_samples, self.n_features_ = X.shape

        is_classification = isinstance(self, ClassifierMixin)

        y = np.atleast_1d(y)
        if y.ndim == 1:
            y = y[:, np.newaxis]

        self.classes_ = []
        self.n_classes_ = []
        self.n_outputs_ = y.shape[1]

        if is_classification:
            y = np.copy(y)

            for k in xrange(self.n_outputs_):
                unique = np.unique(y[:, k])
                self.classes_.append(unique)
                self.n_classes_.append(unique.shape[0])
                y[:, k] = np.searchsorted(unique, y[:, k])

        else:
            self.classes_ = [None] * self.n_outputs_
            self.n_classes_ = [1] * self.n_outputs_

        if getattr(y, "dtype", None) != DOUBLE or not y.flags.contiguous:
            y = np.ascontiguousarray(y, dtype=DOUBLE)

        if is_classification:
            criterion = CLASSIFICATION[self.criterion](self.n_outputs_,
                                                       self.n_classes_)
        else:
            criterion = REGRESSION[self.criterion](self.n_outputs_)

        # Check parameters
        max_depth = np.inf if self.max_depth is None else self.max_depth

        if isinstance(self.max_features, basestring):
            if self.max_features == "auto":
                if is_classification:
                    max_features = max(1, int(np.sqrt(self.n_features_)))
                else:
                    max_features = self.n_features_
            elif self.max_features == "sqrt":
                max_features = max(1, int(np.sqrt(self.n_features_)))
            elif self.max_features == "log2":
                max_features = max(1, int(np.log2(self.n_features_)))
            else:
                raise ValueError(
                    'Invalid value for max_features. Allowed string '
                    'values are "auto", "sqrt" or "log2".')
        elif self.max_features is None:
            max_features = self.n_features_
        else:
            max_features = self.max_features

        if len(y) != n_samples:
            raise ValueError("Number of labels=%d does not match "
                             "number of samples=%d" % (len(y), n_samples))
        if self.min_samples_split <= 0:
            raise ValueError("min_samples_split must be greater than zero.")
        if self.min_samples_leaf <= 0:
            raise ValueError("min_samples_leaf must be greater than zero.")
        if max_depth <= 0:
            raise ValueError("max_depth must be greater than zero. ")
        if self.min_density < 0.0 or self.min_density > 1.0:
            raise ValueError("min_density must be in [0, 1]")
        if not (0 < max_features <= self.n_features_):
            raise ValueError("max_features must be in (0, n_features]")
        if sample_mask is not None and len(sample_mask) != n_samples:
            raise ValueError("Length of sample_mask=%d does not match "
                             "number of samples=%d" % (len(sample_mask),
                                                       n_samples))
        if X_argsorted is not None and len(X_argsorted) != n_samples:
            raise ValueError("Length of X_argsorted=%d does not match "
                             "number of samples=%d" % (len(X_argsorted),
                                                       n_samples))

        # Build tree
        self.tree_ = _tree.Tree(self.n_features_, self.n_classes_,
                                self.n_outputs_, criterion, max_depth,
                                self.min_samples_split, self.min_samples_leaf,
                                self.min_density, max_features,
                                self.find_split_, self.random_state)

        self.tree_.build(X, y, sample_mask=sample_mask,
                         X_argsorted=X_argsorted)

        if self.compute_importances:
            self.feature_importances_ = \
                self.tree_.compute_feature_importances()

        return self

    def predict(self, X):
        """Predict class or regression target for X.

        For a classification model, the predicted class for each sample in X is
        returned. For a regression model, the predicted value based on X is
        returned.

        Parameters
        ----------
        X : array-like of shape = [n_samples, n_features]
            The input samples.

        Returns
        -------
        y : array of shape = [n_samples] or [n_samples, n_outputs]
            The predicted classes, or the predict values.
        """
        if getattr(X, "dtype", None) != DTYPE or X.ndim != 2:
            X = array2d(X, dtype=DTYPE, order="F")

        n_samples, n_features = X.shape

        if self.tree_ is None:
            raise Exception("Tree not initialized. Perform a fit first")

        if self.n_features_ != n_features:
            raise ValueError("Number of features of the model must "
                             " match the input. Model n_features is %s and "
                             " input n_features is %s "
                             % (self.n_features_, n_features))

        P = self.tree_.predict(X)

        if isinstance(self, ClassifierMixin):
            predictions = np.zeros((n_samples, self.n_outputs_))

            for k in xrange(self.n_outputs_):
                predictions[:, k] = self.classes_[k].take(np.argmax(P[:, k],
                                                                    axis=1),
                                                          axis=0)
        else:
            predictions = P[:, :, 0]

        if self.n_outputs_ == 1:
            predictions = predictions.reshape((n_samples, ))

        return predictions


class DecisionTreeClassifier(BaseDecisionTree, ClassifierMixin):
    """A decision tree classifier.

    Parameters
    ----------
    criterion : string, optional (default="gini")
        The function to measure the quality of a split. Supported criteria are
        "gini" for the Gini impurity and "entropy" for the information gain.

    max_depth : integer or None, optional (default=None)
        The maximum depth of the tree. If None, then nodes are expanded until
        all leaves are pure or until all leaves contain less than
        min_samples_split samples.

    min_samples_split : integer, optional (default=1)
        The minimum number of samples required to split an internal node.

    min_samples_leaf : integer, optional (default=1)
        The minimum number of samples required to be at a leaf node.

    min_density : float, optional (default=0.1)
        This parameter controls a trade-off in an optimization heuristic. It
        controls the minimum density of the `sample_mask` (i.e. the
        fraction of samples in the mask). If the density falls below this
        threshold the mask is recomputed and the input data is packed
        which results in data copying.  If `min_density` equals to one,
        the partitions are always represented as copies of the original
        data. Otherwise, partitions are represented as bit masks (aka
        sample masks).

    max_features : int, string or None, optional (default=None)
        The number of features to consider when looking for the best split.
        If "auto", then `max_features=sqrt(n_features)` on classification
        tasks and `max_features=n_features` on regression problems. If "sqrt",
        then `max_features=sqrt(n_features)`. If "log2", then
        `max_features=log2(n_features)`. If None, then
        `max_features=n_features`.

    compute_importances : boolean, optional (default=True)
        Whether feature importances are computed and stored into the
        ``feature_importances_`` attribute when calling fit.

    random_state : int, RandomState instance or None, optional (default=None)
        If int, random_state is the seed used by the random number generator;
        If RandomState instance, random_state is the random number generator;
        If None, the random number generator is the RandomState instance used
        by `np.random`.

    Attributes
    ----------
    `tree_` : Tree object
        The underlying Tree object.

    `feature_importances_` : array of shape = [n_features]
        The feature importances (the higher, the more important the feature).
        The importance I(f) of a feature f is computed as the (normalized)
        total reduction of error brought by that feature. It is also known as
        the Gini importance [4]_.


    See also
    --------
    DecisionTreeRegressor

    References
    ----------

    .. [1] http://en.wikipedia.org/wiki/Decision_tree_learning

    .. [2] L. Breiman, J. Friedman, R. Olshen, and C. Stone, "Classification
           and Regression Trees", Wadsworth, Belmont, CA, 1984.

    .. [3] T. Hastie, R. Tibshirani and J. Friedman. "Elements of Statistical
           Learning", Springer, 2009.

    .. [4] L. Breiman, and A. Cutler, "Random Forests",
           http://www.stat.berkeley.edu/~breiman/RandomForests/cc_home.htm

    Examples
    --------
    >>> from sklearn.datasets import load_iris
    >>> from sklearn.cross_validation import cross_val_score
    >>> from sklearn.tree import DecisionTreeClassifier

    >>> clf = DecisionTreeClassifier(random_state=0)
    >>> iris = load_iris()

    >>> cross_val_score(clf, iris.data, iris.target, cv=10)
    ...                             # doctest: +SKIP
    ...
    array([ 1.     ,  0.93...,  0.86...,  0.93...,  0.93...,
            0.93...,  0.93...,  1.     ,  0.93...,  1.      ])
    """
    def __init__(self, criterion="gini",
                       max_depth=None,
                       min_samples_split=1,
                       min_samples_leaf=1,
                       min_density=0.1,
                       max_features=None,
                       compute_importances=False,
                       random_state=None):
        super(DecisionTreeClassifier, self).__init__(criterion,
                                                     max_depth,
                                                     min_samples_split,
                                                     min_samples_leaf,
                                                     min_density,
                                                     max_features,
                                                     compute_importances,
                                                     random_state)

    def predict_proba(self, X):
        """Predict class probabilities of the input samples X.

        Parameters
        ----------
        X : array-like of shape = [n_samples, n_features]
            The input samples.

        Returns
        -------
        p : array of shape = [n_samples, n_classes], or a list of n_outputs
            such arrays if n_outputs > 1.
            The class probabilities of the input samples. Classes are ordered
            by arithmetical order.
        """
        if getattr(X, "dtype", None) != DTYPE or X.ndim != 2:
            X = array2d(X, dtype=DTYPE, order="F")

        n_samples, n_features = X.shape

        if self.tree_ is None:
            raise Exception("Tree not initialized. Perform a fit first.")

        if self.n_features_ != n_features:
            raise ValueError("Number of features of the model must "
                             " match the input. Model n_features is %s and "
                             " input n_features is %s "
                             % (self.n_features_, n_features))

        proba = []
        P = self.tree_.predict(X)

        for k in xrange(self.n_outputs_):
            P_k = P[:, k, :self.n_classes_[k]]
            normalizer = P_k.sum(axis=1)[:, np.newaxis]
            normalizer[normalizer == 0.0] = 1.0
            P_k /= normalizer
            proba.append(P_k)

        if self.n_outputs_ == 1:
            return proba[0]

        else:
            return proba

    def predict_log_proba(self, X):
        """Predict class log-probabilities of the input samples X.

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


class DecisionTreeRegressor(BaseDecisionTree, RegressorMixin):
    """A tree regressor.

    Parameters
    ----------
    criterion : string, optional (default="mse")
        The function to measure the quality of a split. The only supported
        criterion is "mse" for the mean squared error.

    max_depth : integer or None, optional (default=None)
        The maximum depth of the tree. If None, then nodes are expanded until
        all leaves are pure or until all leaves contain less than
        min_samples_split samples.

    min_samples_split : integer, optional (default=1)
        The minimum number of samples required to split an internal node.

    min_samples_leaf : integer, optional (default=1)
        The minimum number of samples required to be at a leaf node.

    min_density : float, optional (default=0.1)
        This parameter controls a trade-off in an optimization heuristic. It
        controls the minimum density of the `sample_mask` (i.e. the
        fraction of samples in the mask). If the density falls below this
        threshold the mask is recomputed and the input data is packed
        which results in data copying.  If `min_density` equals to one,
        the partitions are always represented as copies of the original
        data. Otherwise, partitions are represented as bit masks (aka
        sample masks).

    max_features : int, string or None, optional (default=None)
        The number of features to consider when looking for the best split.
        If "auto", then `max_features=sqrt(n_features)` on classification
        tasks and `max_features=n_features` on regression problems. If "sqrt",
        then `max_features=sqrt(n_features)`. If "log2", then
        `max_features=log2(n_features)`. If None, then
        `max_features=n_features`.

    compute_importances : boolean, optional (default=True)
        Whether feature importances are computed and stored into the
        ``feature_importances_`` attribute when calling fit.

    random_state : int, RandomState instance or None, optional (default=None)
        If int, random_state is the seed used by the random number generator;
        If RandomState instance, random_state is the random number generator;
        If None, the random number generator is the RandomState instance used
        by `np.random`.

    Attributes
    ----------
    `tree_` : Tree object
        The underlying Tree object.

    `feature_importances_` : array of shape = [n_features]
        The feature importances (the higher, the more important the feature).
        The importance I(f) of a feature f is computed as the (normalized)
        total reduction of error brought by that feature. It is also known as
        the Gini importance [4]_.


    See also
    --------
    DecisionTreeClassifier

    References
    ----------

    .. [1] http://en.wikipedia.org/wiki/Decision_tree_learning

    .. [2] L. Breiman, J. Friedman, R. Olshen, and C. Stone, "Classification
           and Regression Trees", Wadsworth, Belmont, CA, 1984.

    .. [3] T. Hastie, R. Tibshirani and J. Friedman. "Elements of Statistical
           Learning", Springer, 2009.

    .. [4] L. Breiman, and A. Cutler, "Random Forests",
           http://www.stat.berkeley.edu/~breiman/RandomForests/cc_home.htm

    Examples
    --------
    >>> from sklearn.datasets import load_boston
    >>> from sklearn.cross_validation import cross_val_score
    >>> from sklearn.tree import DecisionTreeRegressor

    >>> boston = load_boston()
    >>> regressor = DecisionTreeRegressor(random_state=0)

    R2 scores (a.k.a. coefficient of determination) over 10-folds CV:

    >>> cross_val_score(regressor, boston.data, boston.target, cv=10)
    ...                    # doctest: +SKIP
    ...
    array([ 0.61..., 0.57..., -0.34..., 0.41..., 0.75...,
            0.07..., 0.29..., 0.33..., -1.42..., -1.77...])
    """
    def __init__(self, criterion="mse",
                       max_depth=None,
                       min_samples_split=1,
                       min_samples_leaf=1,
                       min_density=0.1,
                       max_features=None,
                       compute_importances=False,
                       random_state=None):
        super(DecisionTreeRegressor, self).__init__(criterion,
                                                    max_depth,
                                                    min_samples_split,
                                                    min_samples_leaf,
                                                    min_density,
                                                    max_features,
                                                    compute_importances,
                                                    random_state)


class ExtraTreeClassifier(DecisionTreeClassifier):
    """An extremely randomized tree classifier.

    Extra-trees differ from classic decision trees in the way they are built.
    When looking for the best split to separate the samples of a node into two
    groups, random splits are drawn for each of the `max_features` randomly
    selected features and the best split among those is chosen. When
    `max_features` is set 1, this amounts to building a totally random
    decision tree.

    Warning: Extra-trees should only be used within ensemble methods.

    See also
    --------
    ExtraTreeRegressor, ExtraTreesClassifier, ExtraTreesRegressor

    References
    ----------

    .. [1] P. Geurts, D. Ernst., and L. Wehenkel, "Extremely randomized trees",
           Machine Learning, 63(1), 3-42, 2006.
    """
    def __init__(self, criterion="gini",
                       max_depth=None,
                       min_samples_split=1,
                       min_samples_leaf=1,
                       min_density=0.1,
                       max_features="auto",
                       compute_importances=False,
                       random_state=None):
        super(ExtraTreeClassifier, self).__init__(criterion,
                                                  max_depth,
                                                  min_samples_split,
                                                  min_samples_leaf,
                                                  min_density,
                                                  max_features,
                                                  compute_importances,
                                                  random_state)

        self.find_split_ = _tree.TREE_SPLIT_RANDOM


class ExtraTreeRegressor(DecisionTreeRegressor):
    """An extremely randomized tree regressor.

    Extra-trees differ from classic decision trees in the way they are built.
    When looking for the best split to separate the samples of a node into two
    groups, random splits are drawn for each of the `max_features` randomly
    selected features and the best split among those is chosen. When
    `max_features` is set 1, this amounts to building a totally random
    decision tree.

    Warning: Extra-trees should only be used within ensemble methods.

    See also
    --------
    ExtraTreeClassifier : A classifier base on extremely randomized trees
    sklearn.ensemble.ExtraTreesClassifier : An ensemble of extra-trees for
        classification
    sklearn.ensemble.ExtraTreesRegressor : An ensemble of extra-trees for
        regression

    References
    ----------

    .. [1] P. Geurts, D. Ernst., and L. Wehenkel, "Extremely randomized trees",
           Machine Learning, 63(1), 3-42, 2006.
    """
    def __init__(self, criterion="mse",
                       max_depth=None,
                       min_samples_split=1,
                       min_samples_leaf=1,
                       min_density=0.1,
                       max_features="auto",
                       compute_importances=False,
                       random_state=None):
        super(ExtraTreeRegressor, self).__init__(criterion,
                                                 max_depth,
                                                 min_samples_split,
                                                 min_samples_leaf,
                                                 min_density,
                                                 max_features,
                                                 compute_importances,
                                                 random_state)

        self.find_split_ = _tree.TREE_SPLIT_RANDOM
