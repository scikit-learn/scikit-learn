# Adapted from MILK: Machine Learning Toolkit
# Copyright (C) 2008-2011, Luis Pedro Coelho <luis@luispedro.org>
# License: MIT. See COPYING.MIT file in the milk distribution
#
# Authors: Brian Holt, Peter Prettenhofer, Satrajit Ghosh
#
# License: BSD Style.

from __future__ import division
from ..utils import check_random_state
import numpy as np
from ..base import BaseEstimator, ClassifierMixin, RegressorMixin
import _tree

__all__ = [
    'DecisionTreeClassifier',
    'DecisionTreeRegressor',
]

CLASSIFICATION = {
    'gini': _tree.Gini,
    'entropy': _tree.Entropy,
    #'miss': _tree.eval_miss,
}

REGRESSION = {
    'mse': _tree.MSE,
}

GRAPHVIZ_TREE_TEMPLATE = """\
%(tree)s [label="%(tree_gv)s"] ;
[label="%(tree_left_gv)s"] ;
[label="%(tree_right_gv)s"] ;
%(tree)s -> %(tree_left)s ;
%(tree)s -> %(tree_right)s ;

"""


class Leaf(object):
    """A class to store leaf values in the tree.

    Parameters
    ----------

    value : array-like, shape = [n_features] OR 1
        For classification it is a histogram of target values
        For regression is it the mean for the region

    See also
    --------

    Node
    """

    def __init__(self, value):
        self.value = value

    def _graphviz(self):
        """Print the leaf for graph visualisation."""
        return 'Leaf(%s)' % (self.value)


class Node(object):
    """A class to store node information in the tree.

    Parameters
    ----------

    dimension : integer
        The dimension used to split on
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

    See also
    --------

    Leaf
    """

    def __init__(self, dimension, threshold, error, samples, value,
                 left, right):
        self.dimension = dimension
        self.threshold = threshold
        self.error = error
        self.samples = samples
        self.value = value
        self.left = left
        self.right = right

    def _graphviz(self):
        """Print the node for graph visualisation."""

        return "x[%s] < %s \\n error = %s \\n samples = %s \\n v = %s" \
               % (self.dimension, self.threshold,\
                  self.error, self.samples, self.value)


def _build_tree(is_classification, X, y, criterion,
               max_depth, min_split, max_features, n_classes, random_state):
    """Build a tree by recursively partitioning the data."""

    n_samples, n_features = X.shape
    if len(y) != len(X):
        raise ValueError("Number of labels=%d does not match "
                          "number of features=%d"
                         % (len(y), len(X)))
    y = np.array(y, dtype=np.float64, order="c")

    feature_mask = np.ones(n_features, dtype=np.bool)
    sample_dims = np.arange(n_features)
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

    X = X[:, feature_mask]

    # make data fortran layout
    if not X.flags["F_CONTIGUOUS"]:
        X = np.array(X, order="F")

    if min_split <= 0:
        raise ValueError("min_split must be greater than zero. "
                         "min_split is %s." % min_split)
    if max_depth <= 0:
        raise ValueError("max_depth must be greater than zero. "
                         "max_depth is %s." % max_depth)

    def recursive_partition(X, y, depth):
        is_split_valid = True

        if depth >= max_depth:
            is_split_valid = False

        dim, threshold, error, init_error = _tree._find_best_split(
            X, y, criterion)

        if dim != -1:
            split = X[:, dim] < threshold
            if len(X[split]) < min_split or len(X[~split]) < min_split:
                is_split_valid = False
        else:
            is_split_valid = False

        if is_classification:
            a = np.zeros((n_classes,))
            t = y.max() + 1
            a[:t] = np.bincount(y.astype(np.int))
        else:
            a = np.mean(y)

        if is_split_valid == False:
            return Leaf(a)

        left_partition = recursive_partition(X[split], y[split], depth + 1)
        right_partition = recursive_partition(X[~split], y[~split], depth + 1)

        return Node(dimension=sample_dims[dim], threshold=threshold,
                    error=init_error, samples=len(y), value=a,
                    left=left_partition, right=right_partition)

    return recursive_partition(X, y, 0)


def _apply_tree(tree, X):
    """Applies the decision tree to X."""

    if type(tree) is Leaf:
        return tree.value
    if X[tree.dimension] < tree.threshold:
        return _apply_tree(tree.left, X)
    return _apply_tree(tree.right, X)


def _graphviz(tree):
    """Print decision tree in .dot format."""

    if type(tree) is Leaf:
        return ""
    s = GRAPHVIZ_TREE_TEMPLATE % {
        "tree": tree,
        "tree_gv": tree._graphviz(),
        "tree_left": tree.left,
        "tree_left_gv": tree.left._graphviz(),
        "tree_right": tree.right,
        "tree_right_gv": tree.right._graphviz(),
    }
    return s + _graphviz(tree.left) + _graphviz(tree.right)


class BaseDecisionTree(BaseEstimator):
    """Should not be used directly, use derived classes instead."""

    _tree_types = ['classification', 'regression']
    _classification_subtypes = ['binary', 'multiclass']

    def __init__(self, n_classes, impl, criterion, max_depth,
                 min_split, max_features, random_state):

        if not impl in self._tree_types:
            raise ValueError("impl should be one of %s, %s was given"
                             % (self._tree_types, impl))

        self.type = impl
        self.n_classes = n_classes
        self.classification_subtype = None
        self.criterion = criterion
        self.max_depth = max_depth
        self.min_split = min_split
        self.max_features = max_features
        self.random_state = check_random_state(random_state)

        self.n_features = None
        self.tree = None

    def export_to_graphviz(self, filename="tree.dot"):
        """Export the tree in ".dot" format.

        Once exported, you can render to PostScript using, for example,
        $ dot -Tps tree.dot -o tree.ps

        Parameters
        ----------
        filename : string
            The name of the file to write to.

        """
        if self.tree is None:
            raise Exception('Tree not initialized. Perform a fit first')

        with open(filename, 'w') as f:
            f.write("digraph Tree {\n")
            f.write(_graphviz(self.tree))
            f.write("\n}\n")

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
        X = np.asanyarray(X, dtype=np.float64, order='C')
        n_samples, self.n_features = X.shape

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
            label_counts_left = np.zeros((self.n_classes,), dtype=np.int32)
            label_counts_right = np.zeros((self.n_classes,), dtype=np.int32)
            criterion = criterion_class(self.n_classes, label_counts_left,
                                        label_counts_right)

            self.tree = _build_tree(True, X, y, criterion, self.max_depth,
                                    self.min_split, self.max_features,
                                    self.n_classes, self.random_state)
        else:  # regression
            y = np.asanyarray(y, dtype=np.float64, order='C')

            criterion_class = REGRESSION[self.criterion]
            y_temp = np.zeros((n_samples,), dtype=np.float64)
            criterion = criterion_class(y_temp)
            self.tree = _build_tree(False, X, y, criterion, self.max_depth,
                                    self.min_split, self.max_features,
                                    None, self.random_state)
        return self

    def predict(self, X):
        """Predict class or regression target for a test vector X.

        For a classification model, the predicted class for each
        sample in X is returned.  For a regression model, the predicted
        value based on X is returned.

        Parameters
        ----------
        X : array-like, shape = [n_samples, n_features]

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
    >>> import numpy as np
    >>> from sklearn.datasets import load_iris
    >>> from sklearn.cross_val import StratifiedKFold
    >>> from sklearn.tree import DecisionTreeClassifier
    >>> data = load_iris()
    >>> skf = StratifiedKFold(data.target, 10)
    >>> for train_index, test_index in skf:
    ...     clf = DecisionTreeClassifier(random_state=0)
    ...     clf = clf.fit(data.data[train_index], data.target[train_index])
    ...     print np.mean(clf.predict(data.data[test_index]) ==
    ...          data.target[test_index]) #doctest: +ELLIPSIS
    ...
    1.0
    0.933333333333
    0.866666666667
    0.933333333333
    0.933333333333
    0.933333333333
    0.933333333333
    1.0
    0.933333333333
    1.0
    """

    def __init__(self, n_classes=None, criterion='gini', max_depth=10,
                  min_split=1, max_features=None, random_state=None):
        BaseDecisionTree.__init__(self, n_classes, 'classification',
                                  criterion, max_depth, min_split,
                                  max_features, random_state)

    def predict_proba(self, X):
        """Predict class probabilities on a test vector X.

        Parameters
        ----------
        X : array-like, shape = [n_samples, n_features]

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
        """Predict class log probabilities on a test vector X

        Parameters
        ----------
        X : array-like, shape = [n_samples, n_features]

        Returns
        -------
        P : array-like, shape = [n_samples, n_classes]
            Returns the log-probabilities of the sample for each class in
            the model, where classes are ordered by arithmetical
            order.

        """
        return np.log(self.predict_proba(X))


class DecisionTreeRegressor(BaseDecisionTree, RegressorMixin):
    """A tree regressor

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
    >>> import numpy as np
    >>> from sklearn.datasets import load_boston
    >>> from sklearn.cross_val import KFold
    >>> from sklearn.tree import DecisionTreeRegressor
    >>> data = load_boston()
    >>> kf = KFold(len(data.target), 10)
    >>> for train_index, test_index in kf:
    ...     clf = DecisionTreeRegressor(random_state=0)
    ...     clf = clf.fit(data.data[train_index], data.target[train_index])
    ...     print np.mean(np.power(clf.predict(data.data[test_index]) -
    ...          data.target[test_index], 2)) #doctest: +ELLIPSIS
    ...
    12.9450419508
    11.6925868725
    12.8940290384
    59.7824284864
    19.3208876032
    64.0553094769
    15.1038466202
    92.2104637727
    54.4061950617
    50.6928172067
    """

    def __init__(self, criterion='mse', max_depth=10,
                  min_split=1, max_features=None, random_state=None):
        BaseDecisionTree.__init__(self, None, 'regression',
                                  criterion, max_depth, min_split,
                                  max_features, random_state)
