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

import _tree

__all__ = [
    'DecisionTreeClassifier',
    'DecisionTreeRegressor',
]

DTYPE = _tree.DTYPE

CLASSIFICATION = {
    'gini': _tree.Gini,
    'entropy': _tree.Entropy,
    #'miss': _tree.Miss,
}

REGRESSION = {
    'mse': _tree.MSE,
}

GRAPHVIZ_TREE_TEMPLATE = """\
%(current)s [label="%(current_gv)s"] ;
%(left_child)s [label="%(left_child_gv)s"] ;
%(right_child)s [label="%(right_child_gv)s"] ;
%(current)s -> %(left_child)s ;
%(current)s -> %(right_child)s ;
"""


class GraphvizExporter(object):
    """Export a Decision Tree in ".dot" format.

    Once exported, you can render to PostScript using, for example,
    $ dot -Tps tree.dot -o tree.ps

    or to PNG using
    $ dot -Tpng tree.dot -o tree.png
    """

    def __init__(self, out=open("tree.dot", 'w'), feature_names=None):
        """Export a Decision Tree to GraphViz format.

        Generates GraphViz representation of the decision tree. The output
        is written to `out`.

        Parameters
        ----------
        out : file object, optional
            Handle to the output file.

        feature_names : list of strings, optional
            Names of each of the features.
        """
        self.out = out
        self.feature_names = feature_names

    def make_node_repr(self, node):
        if node.is_leaf:
            return "error = %s \\n samples = %s \\n v = %s" \
                % (node.error, node.samples, node.value)
        else:
            feature = "X[%s]" % node.feature
            if self.feature_names is not None:
                feature = self.feature_names[node.feature]

            return "%s < %s \\n error = %s \\n samples = %s \\n v = %s" \
                   % (feature, node.threshold,\
                      node.error, node.samples, node.value)

    def export(self, node):
        """Print the node for graph visualisation.

        Returns
        -------
        self.out : file object
            The file object to which the tree was exported.  The user is
            expected to `close()` this object when done with it.
        """

        def recurse(node):
            current_repr = self.make_node_repr(node)
            left_repr = self.make_node_repr(node.left)
            right_repr = self.make_node_repr(node.right)
            node_data = {
                "current": node.id,
                "current_gv": current_repr,
                "left_child": node.left.id,
                "left_child_gv": left_repr,
                "right_child": node.right.id,
                "right_child_gv": right_repr,
                }
            self.out.write(GRAPHVIZ_TREE_TEMPLATE % node_data)

            if not node.left.is_leaf:
                recurse(node.left)
            if not node.right.is_leaf:
                recurse(node.right)

        self.out.write("digraph Tree {\n")
        recurse(node)
        self.out.write("}")

        return self.out


class Node(object):
    """A class to store node information in the tree.

    Parameters
    ----------

    feature : integer
        The feature used to split on

    threshold : float
        The threshold value to split on

    error : float
        The error in the node.  This could be the impurity (calculated using
        an entropy measure for classification) or the residual regression
        error (calculated using an estimator)

    samples : integer
        Number of samples present at this node

    value : array-like, shape = [n_features] OR 1
        For classification it is a histogram of target values
        For regression is it the mean for the region

    left : Node
        The left child node

    right : Node
        The right child node
    """
    class_counter = 0

    def __init__(self, feature=None, threshold=None, error=None, samples=None,
                 value=None, left=None, right=None):
        self.feature = feature
        self.threshold = threshold
        self.error = error
        self.samples = samples
        self.value = value
        self.left = left
        self.right = right

        if left == None and right == None:
            self.is_leaf = True
        else:
            self.is_leaf = False

        self.id = Node.class_counter
        Node.class_counter += 1


def _build_tree(is_classification, X, y, criterion, max_depth, min_split,
                max_features, n_classes, random_state, min_density,
                sample_mask=None, X_argsorted=None):
    """Build a tree by recursively partitioning the data."""
    assert X.shape[0] == y.shape[0]
    assert max_depth > 0
    assert min_split > 0

    # make data fortran layout
    if not np.isfortran(X):
        X = np.array(X, order="F")

    y = np.array(y, dtype=DTYPE, order="C")

    if X_argsorted is None:
        X_argsorted = np.asfortranarray(
            np.argsort(X.T, axis=1).astype(np.int32).T)

    if sample_mask is None:
        sample_mask = np.ones((X.shape[0],), dtype=np.bool)

    # get num samples from sample_mask instead of X
    n_samples = sample_mask.sum()
    n_features = X.shape[1]

    feature_mask = np.ones((n_features,), dtype=np.bool, order="c")
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

    Node.class_counter = 0

    def recursive_partition(X, X_argsorted, y, sample_mask, depth):
        is_split_valid = True
        n_samples = sample_mask.sum()
        if depth >= max_depth or n_samples < min_split:
            is_split_valid = False
        else:
            feature, threshold, init_error = _tree._find_best_split(
                X, y, X_argsorted, sample_mask, feature_mask,
                criterion, n_samples)

            if feature == -1:
                is_split_valid = False

        current_y = y[sample_mask]
        if is_classification:
            a = np.zeros((n_classes,))
            t = current_y.max() + 1
            a[:t] = np.bincount(current_y.astype(np.int))
        else:
            a = np.mean(current_y)

        if not is_split_valid:
            return Node(error=0.0, samples=n_samples, value=a)
        else:
            if n_samples / X.shape[0] <= min_density:
                # sample_mask too sparse - pack X and X_argsorted
                X = X[sample_mask]
                X_argsorted = np.asfortranarray(
                    np.argsort(X.T, axis=1).astype(np.int32).T)
                y = current_y
                sample_mask = np.ones((X.shape[0],), dtype=np.bool)

            split = X[:, feature] < threshold
            left_partition = recursive_partition(X, X_argsorted, y,
                                                 split & sample_mask,
                                                 depth + 1)
            right_partition = recursive_partition(X, X_argsorted, y,
                                                  ~split & sample_mask,
                                                  depth + 1)

            return Node(feature=feature, threshold=threshold,
                        error=init_error, samples=n_samples, value=a,
                        left=left_partition, right=right_partition)

    return recursive_partition(X, X_argsorted, y, sample_mask, 0)


def _apply_tree(node, x):
    """Applies the decision tree to sample `x`."""
    while node is not None:
        if node.is_leaf:
            return node.value
        if x[node.feature] < node.threshold:
            node = node.left
        else:
            node = node.right


class BaseDecisionTree(BaseEstimator):
    """Should not be used directly, use derived classes instead."""

    _tree_types = ['classification', 'regression']
    _classification_subtypes = ['binary', 'multiclass']

    def __init__(self, n_classes, impl, criterion, max_depth,
                 min_split, max_features, random_state, min_density):

        if not impl in self._tree_types:
            raise ValueError("impl should be one of %s, %s was given"
                             % (self._tree_types, impl))

        self.type = impl
        self.n_classes = n_classes
        self.classification_subtype = None
        self.criterion = criterion

        if min_split <= 0:
            raise ValueError("min_split must be greater than zero.")
        self.min_split = min_split
        if max_depth <= 0:
            raise ValueError("max_depth must be greater than zero. ")
        self.max_depth = max_depth

        self.max_features = max_features
        self.random_state = check_random_state(random_state)

        if min_density < 0.0 or min_density > 1.0:
            raise ValueError("min_density must be in [0, 1]")
        self.min_density = min_density
        self.n_features = None
        self.tree = None

    def export(self, exporter):
        """Export the tree using an exporter.

        Parameters
        ----------
        exporter : class
            Any class that has `export` implemented.

        Example
        -------
        >>> from sklearn.datasets import load_iris
        >>> from sklearn import tree

        >>> clf = tree.DecisionTreeClassifier()
        >>> iris = load_iris()

        >>> clf = clf.fit(iris.data, iris.target)
        >>> import tempfile
        >>> t = tempfile.TemporaryFile()
        >>> exporter = tree.GraphvizExporter(out=t)
        >>> clf.export(exporter)
        >>> t.close()

        """
        if self.tree is None:
            raise Exception('Tree not initialized. Perform a fit first')

        exporter.export(self.tree)

    def fit(self, X, y):
        """Fit the tree with the given training data and parameters.

        Parameters
        ----------
        X : array-like, shape = [n_samples, n_features]
            Training vectors, where n_samples is the number of samples
            and n_features is the number of features.

        y : array-like, shape = [n_samples]
            Target values (integers in classification, real numbers in
            regression)
            For classification, labels must correspond to classes
            0, 1, ..., n_classes-1

        Returns
        -------
        self : object
            Returns self.
        """
        X = np.asanyarray(X, dtype=DTYPE, order='F')
        n_samples, self.n_features = X.shape
        if len(y) != n_samples:
            raise ValueError("Number of labels=%d does not match "
                             "number of features=%d"
                             % (len(y), n_samples))

        sample_mask = np.ones((n_samples,), dtype=np.bool)

        if self.type == 'classification':
            y = np.asanyarray(y, dtype=np.int, order='C')

            y_unique = np.unique(y)
            if tuple(y_unique) == (-1, 1):
                if self.n_classes is None:
                    self.n_classes = 2
                else:
                    if self.n_classes != 2:
                        raise ValueError(
                            "n_classes must equal 2 for binary "
                            "classification: got %d " % self.n_classes)
                self.classification_subtype = "binary"
                y[y == -1] = 0  # normalise target
            elif y.min() >= 0:
                if self.n_classes is None:
                    self.n_classes = y.max() + 1
                else:
                    if self.n_classes < y.max() + 1:
                        raise ValueError("Labels must be in range"
                                         "[0 to %s) " % self.n_classes)
                self.classification_subtype = "multiclass"
            else:
                raise ValueError("Labels must be [-1, 1] for binary and "
                                 "in the range [0 to %s) for multiclass "
                                 "classification " % self.n_classes)

            criterion_class = CLASSIFICATION[self.criterion]
            criterion = criterion_class(self.n_classes)

            self.tree = _build_tree(True, X, y, criterion, self.max_depth,
                                    self.min_split, self.max_features,
                                    self.n_classes, self.random_state,
                                    self.min_density, sample_mask)
        else:  # regression
            y = np.asanyarray(y, dtype=DTYPE, order='C')

            criterion_class = REGRESSION[self.criterion]
            criterion = criterion_class()
            self.tree = _build_tree(False, X, y, criterion, self.max_depth,
                                    self.min_split, self.max_features,
                                    None, self.random_state,
                                    self.min_density, sample_mask)
        return self

    def predict(self, X):
        """Predict class or regression target for a test vector X.

        For a classification model, the predicted class for each
        sample in X is returned.  For a regression model, the predicted
        value based on X is returned.

        Parameters
        ----------
        X : array-like, shape = [n_samples, n_features]
            Test vectors, where n_samples is the number of samples
            and n_features is the number of features.

        Returns
        -------
        predictions : array, shape = [n_samples]
        """

        X = np.atleast_2d(X)
        n_samples, n_features = X.shape

        if self.tree is None:
            raise Exception('Tree not initialized. Perform a fit first')

        if self.n_features != n_features:
            raise ValueError("Number of features of the model must "
                             " match the input. Model n_features is %s and "
                             " input n_features is %s "
                             % (self.n_features, n_features))

        if self.type == 'classification':
            if self.classification_subtype == 'binary':
                predictions = np.zeros(n_samples, dtype=int)
                for idx, sample in enumerate(X):
                    tmp = np.argmax(_apply_tree(self.tree, sample))
                    assert tmp == 0 or tmp == 1
                    predictions[idx] = -1 if tmp == 0 else 1
            elif self.classification_subtype == 'multiclass':
                predictions = np.zeros(n_samples, dtype=int)
                for idx, sample in enumerate(X):
                    predictions[idx] = np.argmax(_apply_tree(self.tree,
                                                             sample))
        else:
            predictions = np.zeros(n_samples, dtype=float)
            for idx, sample in enumerate(X):
                predictions[idx] = _apply_tree(self.tree, sample)

        return predictions


class DecisionTreeClassifier(BaseDecisionTree, ClassifierMixin):
    """A binary or multiclass tree classifier.

    Parameters
    ----------
    n_classes : integer, optional
        number of classes (computed at fit() if not provided)

    criterion : string, optional (default='gini')
        function to measure goodness of split

    max_depth : integer, optional (default=10)
        maximum depth of the tree

    min_split : integer, optional (default=1)
        minimum number of samples required at any leaf node

    max_features : integer, optional (default=None)
        if given, then use a subset (max_features) of features.
        max_features must be in range 0 < max_features <= n_features

    random_state : integer or array_like, optional (default=None)
        seed the random number generator

    min_density : float, optional (default=0.1)
        The minimum density of the sample_mask (i.e. the fraction of samples
        in the mask). If the density falls below this threshold the mask
        is recomputed and the input data is packed which results in data
        copying. If min_density equals one, the partitions are always
        represented as copies of the original data. Otherwise, partitions
        are represented as bit masks (aka sample masks).

    References
    ----------

    http://en.wikipedia.org/wiki/Decision_tree_learning

    L. Breiman, J. Friedman, R. Olshen, and C. Stone. Classification and
    Regression Trees. Wadsworth, Belmont, CA, 1984.

    T. Hastie, R. Tibshirani and J. Friedman.
    Elements of Statistical Learning, Springer, 2009.

    See also
    --------

    DecisionTreeRegressor

    Example
    -------
    >>> from sklearn.datasets import load_iris
    >>> from sklearn.cross_val import cross_val_score
    >>> from sklearn.tree import DecisionTreeClassifier

    >>> clf = DecisionTreeClassifier(random_state=0)
    >>> iris = load_iris()

    >>> cross_val_score(clf, iris.data, iris.target, cv=10)
    ...                             # doctest: +ELLIPSIS +NORMALIZE_WHITESPACE
    ...
    array([ 1.     ,  0.93...,  0.86...,  0.93...,  0.93...,
            0.93...,  0.93...,  1.     ,  0.93...,  1.      ])
    """

    def __init__(self, n_classes=None, criterion='gini', max_depth=10,
                 min_split=1, max_features=None, random_state=None,
                 min_density=0.1):
        BaseDecisionTree.__init__(self, n_classes, 'classification',
                                  criterion, max_depth, min_split,
                                  max_features, random_state,
                                  min_density)

    def predict_proba(self, X):
        """Predict class probabilities on a test vector X.

        Parameters
        ----------
        X : array-like, shape = [n_samples, n_features]
            Test vectors, where n_samples is the number of samples
            and n_features is the number of features.

        Returns
        -------
        P : array-like, shape = [n_samples, n_classes]
            Returns the probability of the sample for each class in
            the model, where classes are ordered by arithmetical
            order.

        """
        X = np.atleast_2d(X)
        n_samples, n_features = X.shape

        if self.tree is None:
            raise Exception('Tree not initialized. Perform a fit first')

        if self.n_features != n_features:
            raise ValueError("Number of features of the model must "
                             " match the input. Model n_features is %s and "
                             " input n_features is %s "
                             % (self.n_features, n_features))

        P = np.zeros((n_samples, self.n_classes))
        for idx, sample in enumerate(X):
            P[idx, :] = _apply_tree(self.tree, sample)
            P[idx, :] /= np.sum(P[idx, :])
        return P

    def predict_log_proba(self, X):
        """Predict class log probabilities on a test vector X.

        Parameters
        ----------
        X : array-like, shape = [n_samples, n_features]
            Test vectors, where n_samples is the number of samples
            and n_features is the number of features.

        Returns
        -------
        P : array-like, shape = [n_samples, n_classes]
            Returns the log-probabilities of the sample for each class in
            the model, where classes are ordered by arithmetical
            order.

        """
        return np.log(self.predict_proba(X))


class DecisionTreeRegressor(BaseDecisionTree, RegressorMixin):
    """A tree regressor.

    Parameters
    ----------
    criterion : string, optional (default='mse')
        function to measure goodness of split

    max_depth : integer, optional (default=10)
        maximum depth of the tree

    min_split : integer, optional (default=1)
        minimum number of samples required at any leaf node

    max_features : integer, optional (default=None)
        if given, then use a subset (max_features) of features.
        max_features must be in range 0 < max_features <= n_features

    random_state : integer or array_like, optional
        seed the random number generator

    min_density : float, optional (default=0.1)
        The minimum density of the sample_mask (i.e. the fraction of samples
        in the mask). If the density falls below this threshold the mask
        is recomputed and the input data is packed which results in data
        copying. If min_density equals one, the partitions are always
        represented as copies of the original data. Otherwise, partitions
        are represented as bit masks (aka sample masks).

    References
    ----------

    http://en.wikipedia.org/wiki/Decision_tree_learning

    L. Breiman, J. Friedman, R. Olshen, and C. Stone. Classification and
    Regression Trees. Wadsworth, Belmont, CA, 1984.

    T. Hastie, R. Tibshirani and J. Friedman.
    Elements of Statistical Learning, Springer, 2009.

    See also
    --------

    DecisionTreeClassifier

    Example
    -------
    >>> from sklearn.datasets import load_boston
    >>> from sklearn.cross_val import cross_val_score
    >>> from sklearn.tree import DecisionTreeRegressor

    >>> boston = load_boston()
    >>> regressor = DecisionTreeRegressor(random_state=0)

    R2 scores (a.k.a. coefficient of determination) over 10-folds CV:

    >>> cross_val_score(regressor, boston.data, boston.target, cv=10)
    ...                    # doctest: +ELLIPSIS +NORMALIZE_WHITESPACE
    ...
    array([ 0.61..., 0.57..., -0.34..., 0.41..., 0.75...,
            0.07..., 0.26..., 0.33..., -1.42..., -1.77...])

    """

    def __init__(self, criterion='mse', max_depth=10,
                 min_split=1, max_features=None, random_state=None,
                 min_density=0.1):
        BaseDecisionTree.__init__(self, None, 'regression',
                                  criterion, max_depth, min_split,
                                  max_features, random_state,
                                  min_density)
