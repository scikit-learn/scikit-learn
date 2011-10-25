# Adapted from MILK: Machine Learning Toolkit
# Copyright (C) 2008-2011, Luis Pedro Coelho <luis@luispedro.org>
# License: MIT. See COPYING.MIT file in the milk distribution
#
# Authors: Brian Holt, Peter Prettenhofer, Satrajit Ghosh
#
# License: BSD Style.

from __future__ import division
import numpy as np

from ..base import BaseEstimator, ClassifierMixin, RegressorMixin
from ..utils import check_random_state

from . import _tree

__all__ = [
    "DecisionTreeClassifier",
    "DecisionTreeRegressor",
]

DTYPE = _tree.DTYPE

CLASSIFICATION = {
    "gini": _tree.Gini,
    "entropy": _tree.Entropy,
}

REGRESSION = {
    "mse": _tree.MSE,
}

GRAPHVIZ_TREE_TEMPLATE = """\
%(current)s [label="%(current_gv)s"] ;
%(left_child)s [label="%(left_child_gv)s"] ;
%(right_child)s [label="%(right_child_gv)s"] ;
%(current)s -> %(left_child)s ;
%(current)s -> %(right_child)s ;
"""


def export_graphviz(decision_tree, out_file=None, feature_names=None):
    """Export a decision tree in DOT format.

    This function generates a GraphViz representation of the decision tree,
    which is then written into `out_file`. Once exported, graphical renderings
    can be generated using, for example,::

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
    >>> out_file = export_graphviz(clf, out_file=tempfile.TemporaryFile())
    >>> out_file.close()
    """
    def node_to_str(node):
        if node.is_leaf:
            return "error = %s\\nsamples = %s\\nvalue = %s" \
                % (node.initial_error, node.samples, node.value)
        else:
            if feature_names is not None:
                feature = feature_names[node.feature]
            else:
                feature = "X[%s]" % node.feature

            return "%s <= %s\\nerror = %s\\nsamples = %s\\nvalue = %s" \
                   % (feature, node.threshold,
                      node.initial_error, node.samples, node.value)

    def recurse(node, count):
        node_data = {
            "current": count,
            "current_gv": node_to_str(node),
            "left_child": 2 * count + 1,
            "left_child_gv": node_to_str(node.left),
            "right_child": 2 * count + 2,
            "right_child_gv": node_to_str(node.right),
        }

        out_file.write(GRAPHVIZ_TREE_TEMPLATE % node_data)

        if not node.left.is_leaf:
            recurse(node.left, 2 * count + 1)
        if not node.right.is_leaf:
            recurse(node.right, 2 * count + 2)

    if out_file is None:
        out_file = open("tree.dot", 'w')
    elif isinstance(out_file, basestring):
        out_file = open(out_file, 'w')

    out_file.write("digraph Tree {\n")
    recurse(decision_tree.tree, 0)
    out_file.write("}")

    return out_file


def _build_tree(is_classification, X, y, criterion, max_depth, min_split,
                max_features, n_classes, random_state, min_density,
                sample_mask=None, X_argsorted=None, store_terminal_region=False):
    """Build a tree by recursively partitioning the data."""
    assert X.shape[0] == y.shape[0]
    assert max_depth > 0
    assert min_split > 0

    # only works if we don't fancy index ever
    if store_terminal_region:
        assert min_density == 0.0

    # make data fortran layout
    if not np.isfortran(X):
        X = np.asfortranarray(X)

    y = np.ascontiguousarray(y, dtype=np.float64)

    if X_argsorted is None:
        X_argsorted = np.asfortranarray(
            np.argsort(X.T, axis=1).astype(np.int32).T)

    if sample_mask is None:
        sample_mask = np.ones((X.shape[0],), dtype=np.bool)

    n_features = X.shape[1]
    
    feature_mask = np.ones((n_features,), dtype=np.bool, order="C")
    if max_features is not None:
        if max_features <= 0 or max_features > n_features:
            raise ValueError("max_features=%d must be in range (0..%d]. "
                             "Did you mean to use None to signal no "
                             "max_features?"
                             % (max_features, n_features))

        permutation = random_state.permutation(n_features)
        sample_dims = np.sort(permutation[-max_features:])
        feature_mask[sample_dims] = False
        feature_mask = np.logical_not(feature_mask)

    feature_mask = feature_mask.astype(np.int32)

    def recursive_partition(X, X_argsorted, y, sample_mask, depth):
        is_split_valid = True
        n_samples = sample_mask.sum()
        if n_samples == 0:
            raise ValueError("Attempting to find a split with an empty sample_mask")
        if depth >= max_depth or n_samples < min_split:
            is_split_valid = False
        else:
            feature, threshold, init_error, best_error = _tree._find_best_split(
                X, y, X_argsorted, sample_mask, feature_mask,
                criterion, n_samples)

            if feature == -1:
                is_split_valid = False

        current_y = y[sample_mask]
        if is_classification:
            value = np.zeros((n_classes,))
            t = current_y.max() + 1
            value[:t] = np.bincount(current_y.astype(np.int))
        else:
            # we need to wrap the mean into an array
            value = np.asanyarray(np.mean(current_y), dtype=np.float64)

        if not is_split_valid:
            # FIXME compute error for leaf
            if store_terminal_region:
                leaf = _tree.Node(-1, 0.0, 0.0, 0.0, n_samples, value, None, None,
                                  sample_mask.nonzero()[0])
            else:
                leaf = _tree.Node(-1, 0.0, 0.0, 0.0, n_samples, value, None, None,
                                  None)
            return leaf
        else:
            if n_samples / X.shape[0] <= min_density:
                # sample_mask too sparse - pack X and X_argsorted
                X = X[sample_mask]
                X_argsorted = np.asfortranarray(
                    np.argsort(X.T, axis=1).astype(np.int32).T)
                y = current_y
                sample_mask = np.ones((X.shape[0],), dtype=np.bool)

            split = X[:, feature] <= threshold

            left_partition = recursive_partition(X, X_argsorted, y,
                                                 split & sample_mask,
                                                 depth + 1)
            right_partition = recursive_partition(X, X_argsorted, y,
                                                  ~split & sample_mask,
                                                  depth + 1)

            return _tree.Node(feature, threshold, init_error, best_error,
                        n_samples, value,
                        left_partition, right_partition, None)

    return recursive_partition(X, X_argsorted, y, sample_mask, 0)


class BaseDecisionTree(BaseEstimator):
    """Warning: This class should not be used directly. Use derived classes
       instead."""

    _tree_types = ["classification", "regression"]

    def __init__(self, impl, criterion, max_depth, min_split,
                 max_features, random_state, min_density):

        if not impl in self._tree_types:
            raise ValueError("impl should be one of %s, %s was given"
                             % (self._tree_types, impl))

        self.type = impl
        self.criterion = criterion
        self.min_split = min_split
        self.max_depth = max_depth
        self.max_features = max_features
        self.random_state = random_state
        self.min_density = min_density

        self.classification_subtype = None
        self.n_features = None
        self.tree = None

    def fit(self, X, y):
        """Build a decision tree from the training set (X, y).

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
        X = np.asanyarray(X, dtype=DTYPE, order='F')
        n_samples, self.n_features = X.shape
        if len(y) != n_samples:
            raise ValueError("Number of labels=%d does not match "
                             "number of features=%d"
                             % (len(y), n_samples))
        random_state = check_random_state(self.random_state)
        if self.min_split <= 0:
            raise ValueError("min_split must be greater than zero.")
        if self.max_depth <= 0:
            raise ValueError("max_depth must be greater than zero. ")
        if self.min_density < 0.0 or self.min_density > 1.0:
            raise ValueError("min_density must be in [0, 1]")

        sample_mask = np.ones((n_samples,), dtype=np.bool)

        is_classification = (self.type == "classification")
        if is_classification:
            y = np.ascontiguousarray(y, dtype=np.int)
            self.classes = np.unique(y)
            self.n_classes = self.classes.shape[0]
            y = np.searchsorted(self.classes, y)

            criterion_class = CLASSIFICATION[self.criterion]
            criterion = criterion_class(self.n_classes)
        else:  # regression
            y = np.ascontiguousarray(y, dtype=np.float64)
            self.n_classes = 1

            criterion_class = REGRESSION[self.criterion]
            criterion = criterion_class()

        self.tree = _build_tree(is_classification, X, y, criterion,
                                self.max_depth, self.min_split,
                                self.max_features, self.n_classes,
                                random_state, self.min_density,
                                sample_mask)
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
        predictions : array of shape = [n_samples]
            The predicted classes, or the predict values.
        """
        X = np.atleast_2d(X)
        X = X.astype(DTYPE)
        n_samples, n_features = X.shape

        if self.tree is None:
            raise Exception('Tree not initialized. Perform a fit first')

        if self.n_features != n_features:
            raise ValueError("Number of features of the model must "
                             " match the input. Model n_features is %s and "
                             " input n_features is %s "
                             % (self.n_features, n_features))

        if self.type == "classification":
            predictions = self.classes.take(np.argmax(
                _tree.apply_tree(self.tree, X, self.n_classes), axis=1), axis=0)
        else:
            predictions = _tree.apply_tree(self.tree, X, self.n_classes)
            predictions = predictions.ravel()

        return predictions


class DecisionTreeClassifier(BaseDecisionTree, ClassifierMixin):
    """A decision tree classifier.

    Parameters
    ----------
    criterion : string, optional (default="gini")
        The function to measure the quality of a split. Supported criteria are
        "gini" for the Gini impurity and "entropy" for the information gain.

    max_depth : integer, optional (default=10)
        The maximum depth of the tree.

    min_split : integer, optional (default=1)
        The minimum number of samples required to split an internal node.

    max_features : integer, optional (default=None)
        If given, then use a subset (max_features) of features.
        max_features must be in range 0 < `max_features` <= `n_features`

    random_state : int, RandomState instance or None, optional (default=None)
        If int, random_state is the seed used by the random number generator;
        If RandomState instance, random_state is the random number generator;
        If None, the random number generator is the RandomState instance used
        by `np.random`.

    min_density : float, optional (default=0.1)
        The minimum density of the `sample_mask` (i.e. the fraction of samples
        in the mask). If the density falls below this threshold the mask is
        recomputed and the input data is packed which results in data copying.
        If `min_density` equals to one, the partitions are always represented
        as copies of the original data. Otherwise, partitions are represented
        as bit masks (aka sample masks).

    References
    ----------
    .. [1] http://en.wikipedia.org/wiki/Decision_tree_learning

    .. [2] L. Breiman, J. Friedman, R. Olshen, and C. Stone, "Classification
           and Regression Trees", Wadsworth, Belmont, CA, 1984.

    .. [3] T. Hastie, R. Tibshirani and J. Friedman. "Elements of Statistical
           Learning", Springer, 2009.

    See also
    --------
    DecisionTreeRegressor

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

    def __init__(self, criterion="gini", max_depth=10, min_split=1,
                 max_features=None, random_state=None, min_density=0.1):
        super(DecisionTreeClassifier, self).__init__(
            "classification", criterion, max_depth, min_split, max_features,
            random_state, min_density)

    def predict_proba(self, X):
        """Predict class probabilities of the input samples X.

        Parameters
        ----------
        X : array-like of shape = [n_samples, n_features]
            The input samples.

        Returns
        -------
        P : array of shape = [n_samples, n_classes]
            The class probabilities of the input samples. Classes are ordered
            by arithmetical order.
        """
        X = np.atleast_2d(X)
        X = X.astype(DTYPE)
        n_samples, n_features = X.shape

        if self.tree is None:
            raise Exception("Tree not initialized. Perform a fit first.")

        if self.n_features != n_features:
            raise ValueError("Number of features of the model must "
                             " match the input. Model n_features is %s and "
                             " input n_features is %s "
                             % (self.n_features, n_features))

        P = _tree.apply_tree(self.tree, X, self.n_classes)
        P /= P.sum(axis=1)[:, np.newaxis]
        return P

    def predict_log_proba(self, X):
        """Predict class log-probabilities of the input samples X.

        Parameters
        ----------
        X : array-like of shape = [n_samples, n_features]
            The input samples.

        Returns
        -------
        P : array of shape = [n_samples, n_classes]
            The class log-probabilities of the input samples. Classes are
            ordered by arithmetical order.
        """
        return np.log(self.predict_proba(X))


class DecisionTreeRegressor(BaseDecisionTree, RegressorMixin):
    """A tree regressor.

    Parameters
    ----------
    criterion : string, optional (default="mse")
        The function to measure the quality of a split. The only supported
        criterion is "mse" for the mean squared error.

    max_depth : integer, optional (default=10)
        The maximum depth of the tree.

    min_split : integer, optional (default=1)
        The minimum number of samples required to split an internal node.

    max_features : integer, optional (default=None)
        If given, then use a subset (max_features) of features.
        max_features must be in range 0 < `max_features` <= `n_features`

    random_state : int, RandomState instance or None, optional (default=None)
        If int, random_state is the seed used by the random number generator;
        If RandomState instance, random_state is the random number generator;
        If None, the random number generator is the RandomState instance used
        by `np.random`.

    min_density : float, optional (default=0.1)
        The minimum density of the `sample_mask` (i.e. the fraction of samples
        in the mask). If the density falls below this threshold the mask is
        recomputed and the input data is packed which results in data copying.
        If `min_density` equals to one, the partitions are always represented
        as copies of the original data. Otherwise, partitions are represented
        as bit masks (aka sample masks).

    References
    ----------
    .. [1] http://en.wikipedia.org/wiki/Decision_tree_learning

    .. [2] L. Breiman, J. Friedman, R. Olshen, and C. Stone, "Classification
           and Regression Trees", Wadsworth, Belmont, CA, 1984.

    .. [3] T. Hastie, R. Tibshirani and J. Friedman. "Elements of Statistical
           Learning", Springer, 2009.

    See also
    --------
    DecisionTreeClassifier

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

    def __init__(self, criterion="mse", max_depth=10, min_split=1,
                 max_features=None, random_state=None, min_density=0.1):
        super(DecisionTreeRegressor, self).__init__(
            "regression", criterion, max_depth, min_split, max_features,
            random_state, min_density)
