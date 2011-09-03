# -*- coding: utf-8 -*-
# Copyright (C) 2008-2011, Luis Pedro Coelho <luis@luispedro.org>
# vim: set ts=4 sts=4 sw=4 expandtab smartindent:
#
# License: MIT. See COPYING.MIT file in the milk distribution
"""
================
Tree Classifier
================

A decision tree classifier

Implements Classification and Regression Trees (Breiman et al. 1984)

"""

from __future__ import division
from ..utils import check_random_state
import numpy as np
from ..base import BaseEstimator, ClassifierMixin, RegressorMixin
import _tree

__all__ = [
    'DecisionTreeClassifier',
    'DecisionTreeRegressor',
    ]

lookup_c = {
       'gini': _tree.Gini,
       'entropy': _tree.Entropy,
       #'miss': _tree.eval_miss,
       }
lookup_r = {
       'mse': _tree.MSE,
       }


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


def _build_tree(is_classification, features, labels, criterion,
               max_depth, min_split, max_features, n_classes, random_state):

    n_samples, n_dims = features.shape
    if len(labels) != len(features):
        raise ValueError("Number of labels=%d does not match "
                          "number of features=%d\n"
                         % (len(labels), len(features)))
    labels = np.array(labels, dtype=np.float64, order="c")

    feature_mask = np.ones(n_dims, dtype=np.bool)
    sample_dims = np.arange(n_dims)
    if max_features is not None:
        if max_features <= 0 or max_features > n_dims:
            raise ValueError("max_features=%d must be in range (0..%d].\n"
                             "Did you mean to use None to signal no "
                             "max_features?"
                             % (max_features, n_dims))

        permutation = random_state.permutation(n_dims)
        sample_dims = np.sort(permutation[-max_features:])
        feature_mask[sample_dims] = False
        feature_mask = np.logical_not(feature_mask)

    features = features[:, feature_mask]

    # make data fortran layout
    if not features.flags["F_CONTIGUOUS"]:
        features = np.array(features, order="F")

    if min_split <= 0:
        raise ValueError("min_split must be greater than zero.\n"
                         "min_split is %s." % min_split)
    if max_depth <= 0:
        raise ValueError("max_depth must be greater than zero.\n"
                         "max_depth is %s." % max_depth)

    def recursive_partition(features, labels, depth):
        is_split_valid = True

        if depth >= max_depth:
            is_split_valid = False

        dim, thresh, error, init_error = _tree._find_best_split(features,
                                                                labels,
                                                                criterion)

        if dim != -1:
            split = features[:, dim] < thresh
            if len(features[split]) < min_split or \
               len(features[~split]) < min_split:
                is_split_valid = False
        else:
            is_split_valid = False

        if is_classification:
            a = np.zeros((n_classes, ))
            t = labels.max() + 1
            a[:t] = np.bincount(labels.astype(np.int))
        else:
            a = np.mean(labels)

        if is_split_valid == False:
            return Leaf(a)

        return Node(dimension=sample_dims[dim],
                    threshold=thresh,
                    error=init_error,
                    samples=len(labels),
                    value=a,
                    left=recursive_partition(features[split],
                                             labels[split], depth + 1),
                    right=recursive_partition(features[~split],
                                              labels[~split], depth + 1))

    return recursive_partition(features, labels, 0)


def _apply_tree(tree, features):
    """Applies the decision tree to features."""

    if type(tree) is Leaf:
        return tree.value
    if features[tree.dimension] < tree.threshold:
        return _apply_tree(tree.left, features)
    return _apply_tree(tree.right, features)


def _graphviz(tree):
    """Print decision tree in .dot format."""

    if type(tree) is Leaf:
        return ""
    s = str(tree) + \
        " [label=" + "\"" + tree._graphviz() + "\"" + "] ;\n"
    s += str(tree.left) + \
        " [label=" + "\"" + tree.left._graphviz() + "\"" + "] ;\n"
    s += str(tree.right) + \
        " [label=" + "\"" + tree.right._graphviz() + "\"" + "] ;\n"

    s += str(tree) + " -> " + str(tree.left) + " ;\n"
    s += str(tree) + " -> " + str(tree.right) + " ;\n"

    return s + _graphviz(tree.left) + _graphviz(tree.right)


class BaseDecisionTree(BaseEstimator):
    '''
    Should not be used directly, use derived classes instead
    '''

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
        """
        Export the tree in .dot format.  Render to PostScript using e.g.
        $ dot -Tps tree.dot -o tree.ps

        Parameters
        ----------
        filename : str
            The name of the file to write to.

        """
        if self.tree is None:
            raise Exception('Tree not initialized. Perform a fit first')

        with open(filename, 'w') as f:
            f.write("digraph Tree {\n")
            f.write(_graphviz(self.tree))
            f.write("\n}\n")

    def fit(self, X, y):
        """
        Fit the tree model according to the given training data and
        parameters.

        Parameters
        ----------
        X : array-like, shape = [n_samples, n_features]
            Training vectors, where n_samples is the number of samples
            and n_features is the number of features.

        y : array-like, shape = [n_samples]
            Target values (integers in classification, real numbers in
            regression)
            For classification, labels must correspond to classes 0,1,...,n_classes-1

        Returns
        -------
        self : object
            Returns self.
        """
        X = np.asanyarray(X, dtype=np.float64, order='C')
        n_samples, self.n_features = X.shape

        if self.type == 'classification':
            y = np.asanyarray(y, dtype=np.int, order='C')

            labels = np.unique(y)
            if tuple(labels) == (-1, 1):
                if self.n_classes is None:
                    self.n_classes = 2
                else:
                    if self.n_classes != 2:
                        raise ValueError("n_classes must equal 2 for binary"
                                         " classification")
                self.classification_subtype = "binary"
                y[y == -1] = 0  # normalise target
            elif labels.min() >= 0:
                if self.n_classes is None:
                    self.n_classes = labels.max() + 1
                else:
                    if self.n_classes < labels.max() + 1:
                        raise ValueError("Labels must be in range"
                                         "[0 to %s) " % self.n_classes)
                self.classification_subtype = "multiclass"
            else:
                raise ValueError("Labels must be [-1, 1] for binary and "
                                 "in the range [0 to %s) for multiclass "
                                 "classification " % self.n_classes)

            criterion_class = lookup_c[self.criterion]
            pm_left = np.zeros((self.n_classes,), dtype=np.int32)
            pm_right = np.zeros((self.n_classes,), dtype=np.int32)
            criterion = criterion_class(self.n_classes, pm_left, pm_right)

            self.tree = _build_tree(True, X, y, criterion,
                                    self.max_depth, self.min_split, self.max_features,
                                    self.n_classes, self.random_state)
        else:  # regression
            y = np.asanyarray(y, dtype=np.float64, order='C')

            criterion_class = lookup_r[self.criterion]
            labels_temp = np.zeros((n_samples,), dtype=np.float64)
            criterion = criterion_class(labels_temp)            
            self.tree = _build_tree(False, X, y, criterion,
                                    self.max_depth, self.min_split, self.max_features,
                                    None, self.random_state)
        return self

    def predict(self, X):
        """
        This function does classification or regression on an array of
        test vectors X.

        For a classification model, the predicted class for each
        sample in X is returned.  For a regression model, the function
        value of X calculated is returned.

        Parameters
        ----------
        X : array-like, shape = [n_samples, n_features]

        Returns
        -------
        C : array, shape = [n_samples]
        """

        X = np.atleast_2d(X)
        n_samples, n_features = X.shape

        if self.tree is None:
            raise Exception('Tree not initialized. Perform a fit first')

        if self.n_features != n_features:
            raise ValueError("Number of features of the model must "
                             " match the input.\n"
                             "Model n_features is %s and "
                             " input n_features is %s "
                             % (self.n_features, n_features))

        if self.type == 'classification':
            if self.classification_subtype == 'binary':
                C = np.zeros(n_samples, dtype=int)
                for idx, sample in enumerate(X):
                    tmp = np.argmax(_apply_tree(self.tree, sample))
                    assert tmp == 0 or tmp == 1
                    C[idx] = -1 if tmp == 0 else 1
            elif self.classification_subtype == 'multiclass':
                C = np.zeros(n_samples, dtype=int)
                for idx, sample in enumerate(X):
                    C[idx] = np.argmax(_apply_tree(self.tree, sample))
        else:
            C = np.zeros(n_samples, dtype=float)
            for idx, sample in enumerate(X):
                C[idx] = _apply_tree(self.tree, sample)

        return C


class DecisionTreeClassifier(BaseDecisionTree, ClassifierMixin):
    """Classify a multi-labeled dataset with a decision tree.

    Parameters
    ----------
    n_classes : integer, optional
        number of classes (computed at fit() if not provided)

    criterion : string, optional (default='gini')
        function to measure goodness of split

    max_depth : integer, optional (default=10)
        maximum depth of the tree

    min_split : integer, optional (default=1)
        minimum size to split on

    max_features : integer, optional (default=None)
        if given, then use a subset (max_features) of features.
        max_features must be in range 0 < max_features <= n_features

    random_state : integer or array_like, optional (default=None)
        seed the random number generator


    Example
    -------
    >>> import numpy as np
    >>> from scikits.learn.datasets import load_iris
    >>> from scikits.learn.cross_val import StratifiedKFold
    >>> from scikits.learn.tree import DecisionTreeClassifier
    >>> data = load_iris()
    >>> skf = StratifiedKFold(data.target, 10)
    >>> for train_index, test_index in skf:
    ...     clf = DecisionTreeClassifier(random_state=0)
    ...     clf = clf.fit(data.data[train_index], data.target[train_index])
    ...     print np.mean(clf.predict(data.data[test_index]) == \
                data.target[test_index]) #doctest: +ELLIPSIS
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
        """
        This function does classification on a test vector X
        given a model with probability information.

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
                             " match the input.\n"
                             "Model n_features is %s and "
                             " input n_features is %s "
                             % (self.n_features, n_features))

        P = np.zeros((n_samples, self.n_classes))
        for idx, sample in enumerate(X):
            P[idx, :] = _apply_tree(self.tree, sample)
            P[idx, :] /= np.sum(P[idx, :])
        return P

    def predict_log_proba(self, X):
        """
        This function does classification on a test vector X

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
    """Perform regression on dataset with a decision tree.

    Parameters
    ----------
    criterion : string, optional (default='mse')
        function to measure goodness of split

    max_depth : integer, optional (default=10)
        maximum depth of the tree

    min_split : integer, optional (default=1)
        minimum size to split on

    max_features : integer, optional (default=None)
        if given, then use a subset (max_features) of features.
        max_features must be in range 0 < max_features <= n_features

    random_state : integer or array_like, optional
        seed the random number generator

    Example
    -------
    >>> import numpy as np
    >>> from scikits.learn.datasets import load_boston
    >>> from scikits.learn.cross_val import KFold
    >>> from scikits.learn.tree import DecisionTreeRegressor
    >>> data = load_boston()
    >>> kf = KFold(len(data.target), 10)
    >>> for train_index, test_index in kf:
    ...     clf = DecisionTreeRegressor(random_state=0)
    ...     clf = clf.fit(data.data[train_index], data.target[train_index])
    ...     print np.mean(np.power(clf.predict(data.data[test_index]) - \
                data.target[test_index], 2)) #doctest: +ELLIPSIS
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
