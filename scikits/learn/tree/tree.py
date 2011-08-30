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
import random

__all__ = [
    'DecisionTreeClassifier',
    'DecisionTreeRegressor',
    ]

lookup_c = \
      {'gini': _tree.Gini,
       'entropy': _tree.Entropy,
       #'miss': _tree.eval_miss,
       }
lookup_r = \
      {'mse': _tree.MSE,
      }


class Leaf(object):
    '''
        v : target value
            Classification: array-like, shape = [n_features]
                Histogram of target values
            Regression:  real number
                Mean for the region
    '''

    def __init__(self, v):
        self.v = v

    def _graphviz(self):
        return 'Leaf(%s)' % (self.v)


class Node(object):
    '''
        value : target value
            Classification: array-like, shape = [n_features]
                Histogram of target values
            Regression:  real number
                Mean for the region
    '''

    def __init__(self, dimension, value, error, samples, v, left, right):
        self.dimension = dimension
        self.value = value
        self.error = error
        self.samples = samples
        self.v = v
        self.left = left
        self.right = right

    def _graphviz(self):
        return "x[%s] < %s \\n error = %s \\n samples = %s \\n v = %s" \
               % (self.dimension, self.value, self.error, self.samples, self.v)


def _build_tree(is_classification, features, labels, criterion,
               max_depth, min_split, F, K, random_state):

    n_samples, n_dims = features.shape
    if len(labels) != len(features):
        raise ValueError("Number of labels does not match "
                          "number of features\n"
                         "num labels is %s and num features is %s "
                         % (len(labels), len(features)))
    labels = np.array(labels, dtype=np.float64, order="c")

    sample_dims = np.arange(n_dims)
    if F is not None:
        if F <= 0:
            raise ValueError("F must be > 0.\n"
                             "Did you mean to use None to signal no F?")
        if F > n_dims:
            raise ValueError("F must be < num dimensions of features.\n"
                             "F is %s, n_dims = %s "
                             % (F, n_dims))

        sample_dims = np.unique(random_state.randint(0, n_dims, F))
        # This ugly construction is required because np.random
        # does not have a sample(population, K) method that will
        # return a subset of K elements from the population.
        # In certain instances, random_state.randint() will return duplicates,
        # meaning that len(sample_dims) < F.  This breaks the contract with
        # the user.
        while len(sample_dims) != F:
            sample_dims = np.unique(random_state.randint(0, n_dims, F))
        features = features[:, sample_dims]

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
            a = np.zeros((K, ))
            t = labels.max() + 1
            a[:t] = np.bincount(labels.astype(np.int))
        else:
            a = np.mean(labels)
                
        if is_split_valid == False:
            return Leaf(a)
                
        return Node(dimension=sample_dims[dim],
                    value=thresh,
                    error=init_error,
                    samples=len(labels),
                    v=a,
                    left=recursive_partition(features[split],
                                             labels[split], depth + 1),
                    right=recursive_partition(features[~split],
                                              labels[~split], depth + 1))

    return recursive_partition(features, labels, 0)


def _apply_tree(tree, features):
    '''
    conf = apply_tree(tree, features)

    Applies the decision tree to a set of features.
    '''
    if type(tree) is Leaf:
        return tree.v
    if features[tree.dimension] < tree.value:
        return _apply_tree(tree.left, features)
    return _apply_tree(tree.right, features)


def _graphviz(tree):
    '''Print decision tree in .dot format
    '''
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

    def __init__(self, K, impl, criterion, max_depth,
                 min_split, F, random_state):

        if not impl in self._tree_types:
            raise ValueError("impl should be one of %s, %s was given"
                             % (self._tree_types, impl))

        self.type = impl
        self.K = K
        self.classification_subtype = None
        self.criterion = criterion
        self.max_depth = max_depth
        self.min_split = min_split
        self.F = F
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
            For classification, labels must correspond to classes 0,1,...,K-1

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
                if self.K is None:
                    self.K = 2
                else:
                    if self.K != 2:
                        raise ValueError("K must equal 2 for binary"
                                         " classification")
                self.classification_subtype = "binary"
                y[y == -1] = 0  # normalise target
            elif labels.min() >= 0:
                if self.K is None:
                    self.K = labels.max() + 1
                else:
                    if self.K != 2:
                        raise ValueError("Labels must be in range"
                                         "[0 to %s) " % self.K)
                self.classification_subtype = "multiclass"
            else:
                raise ValueError("Labels must be [-1, 1] for binary and "
                                 "in the range [0 to %s) for multiclass "
                                 "classification " % self.K)

            criterion_class = lookup_c[self.criterion]
            pm_left = np.zeros((self.K,), dtype=np.int)
            pm_right = np.zeros((self.K,), dtype=np.int)
            criterion = criterion_class(self.K, pm_left, pm_right)

            self.tree = _build_tree(True, X, y, criterion,
                                    self.max_depth, self.min_split, self.F,
                                    self.K, self.random_state)
        else:  # regression
            y = np.asanyarray(y, dtype=np.float64, order='C')
            
            criterion_class = lookup_r[self.criterion]
            labels_temp = np.zeros((n_samples,), dtype=np.float)
            criterion = criterion_class(labels_temp)            
            self.tree = _build_tree(False, X, y, criterion,
                                    self.max_depth, self.min_split, self.F,
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
    K : integer, optional
        number of classes (computed at fit() if not provided)

    criterion : string, optional (default='gini')
        function to measure goodness of split

    max_depth : integer, optional (default=10)
        maximum depth of the tree

    min_split : integer, optional (default=1)
        minimum size to split on

    F : integer, optional (default=None)
        if given, then, choose F features, 0 < F <= n_dims

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
    0.933333333333
    0.933333333333
    0.733333333333
    0.933333333333
    0.933333333333
    0.933333333333
    0.933333333333
    0.933333333333
    0.866666666667
    1.0
    """

    def __init__(self, K=None, criterion='gini', max_depth=10,
                  min_split=1, F=None, random_state=None):
        BaseDecisionTree.__init__(self, K, 'classification', criterion,
                                  max_depth, min_split, F, random_state)

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

        P = np.zeros((n_samples, self.K))
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

    F : integer, optional (default=None)
        if given, then, choose F features, 0 < F <= n_dims

    random_state : integer or array_like, optional
        seed the random number generator

    Example
    -------
    >>> import numpy as np
    >>> from scikits.learn.datasets import load_boston
    >>> from scikits.learn.cross_val import KFold
    >>> from scikits.learn.tree import DecisionTreeRegressor
    >>> data = load_boston()
    >>> kf = KFold(len(data.target), 2)
    >>> for train_index, test_index in kf:
    ...     clf = DecisionTreeRegressor(random_state=0)
    ...     clf = clf.fit(data.data[train_index], data.target[train_index])
    ...     print np.mean(np.power(clf.predict(data.data[test_index]) - \
                data.target[test_index], 2)) #doctest: +ELLIPSIS
    ...
    19.8246021636
    40.4481135756
    """

    def __init__(self, criterion='mse', max_depth=10,
                  min_split=1, F=None, random_state=None):
        BaseDecisionTree.__init__(self, None, 'regression', criterion,
                                  max_depth, min_split, F, random_state)
