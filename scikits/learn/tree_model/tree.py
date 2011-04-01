# -*- coding: utf-8 -*-
# Copyright (C) 2008-2011, Luis Pedro Coelho <luis@luispedro.org>
# vim: set ts=4 sts=4 sw=4 expandtab smartindent:
#
# License: MIT. See COPYING.MIT file in the milk distribution
'''
================
Tree Classifier
================

Decision tree based classifier
------------------------------

'''

from __future__ import division
import numpy as np
from .base import normaliselabels
from .criteria import information_gain
from ..base import BaseEstimator, ClassifierMixin, RegressorMixin

__all__ = [
    'DecisionTreeClassifier',
    'DecisionTreeRegressor',
    ]

class Leaf(object):
    '''
    v : value
    w : weight
    '''
    def __init__(self, v, w):
        self.v = v
        self.w = w
    def __repr__(self):
        return 'Leaf(%s,%s)' % (self.v, self.w)

class Node(object): # This could be replaced by a namedtuple
    def __init__(self, featid, featval, left, right):
        self.featid = featid
        self.featval = featval
        self.left = left
        self.right = right

def _split(features, labels, weights, criterion, subsample, R):
    N,f = features.shape
    if subsample is not None:
        samples = np.array(R.sample(xrange(features.shape[1]), subsample))
        features = features[:, samples]
        f = subsample
    best = None
    best_val = -1.
    for i in xrange(f):
        domain_i = sorted(set(features[:,i]))
        for d in domain_i[1:]:
            cur_split = (features[:,i] < d)
            if weights is not None:
                value = criterion(labels[cur_split], labels[~cur_split], weights[cur_split], weights[~cur_split])
            else:
                value = criterion(labels[cur_split], labels[~cur_split])
            if value > best_val:
                best_val = value
                if subsample is not None:
                    ti = samples[i]
                else:
                    ti = i
                best = ti,d
    return best


def build_tree(features, labels, criterion, min_split=4, subsample=None, R=None, weights=None):
    '''
    tree = build_tree(features, labels, criterion, min_split=4, subsample=None, R=None, weights={all 1s})

    Parameters
    ----------
    features : sequence
        features to use
    labels : sequence
        labels
    criterion : function {labels} x {labels} -> float
        function to measure goodness of split
    min_split : integer
        minimum size to split on
    subsample : integer, optional
        if given, then, at each step, choose subsample features
    R : source of randomness, optional
        See `get_pyrandom`
    weights : sequence, optional
        weight of instance (default: all the same)

    Returns
    -------
    tree : Tree
    '''
    assert len(features) == len(labels), 'build_tree: Nr of labels does not match nr of features'
    features = np.asanyarray(features)
    labels = np.asanyarray(labels)
    if subsample is not None:
        if subsample <= 0:
            raise ValueError('tree_model.build_tree: `subsample` must be > 0.\nDid you mean to use None to signal no subsample?')
        from .base import get_pyrandom
        R = get_pyrandom(R)

    def recursive(features, labels):
        N = float(len(labels))
        if N < min_split:
            return Leaf(labels.sum()/N, N)
        S = _split(features, labels, weights, criterion, subsample, R)
        if S is None:
            return Leaf(labels.sum()/N, N)
        i,thresh = S
        split = features[:,i] < thresh
        return Node(featid=i,
                    featval=thresh,
                    left =recursive(features[ split], labels[ split]),
                    right=recursive(features[~split], labels[~split]))
    return recursive(features, labels)

def apply_tree(tree, features):
    '''
    conf = apply_tree(tree, features)

    Applies the decision tree to a set of features.
    '''
    if type(tree) is Leaf:
        return tree.v
    if features[tree.featid] < tree.featval:
        return apply_tree(tree.left, features)
    return apply_tree(tree.right, features)

def print_tree(tree):
    '''Print decision tree
    '''
    pass

class DecisionTree(BaseEstimator):
    '''A decision tree classifier (currently, implements the greedy ID3 algorithm without any pruning).

    Attributes
    ----------
    
    criterion : function, optional
        criterion to use for tree construction,
        this should be a function that receives two pairs of labels and returns
        a value that indicates how much more information is contained in one
        set of labels over the other.

    min_split : integer, optional
        minimum size to split on (default: 4).
        
    n_subsample : integer, optional
        if given, then, at each step, choose
        
    R : source of randomness, optional
        See `get_pyrandom`
    
    '''
    def __init__(self, criterion=None, min_split=4, n_subsample=None, R=None):
        self.criterion = criterion
        self.min_split = min_split
        self.subsample = n_subsample
        self.R = R
        self.tree = None

    def fit(self, X, y, weights=None):
        self.tree = build_tree(X, y, self.criterion, self.min_split, self.subsample, self.R, weights)
        return self

    def __repr__(self):
        print_tree(self.tree)
        
    def predict(self, X):
        if self.tree is None:
            raise Exception('Tree not initialized. Perform a fit first')
        n_observations,_ = X.shape
        values = np.zeros(n_observations)
        for idx, observation in enumerate(X):
            values[idx] = apply_tree(self.tree, observation)
        return values

class DecisionTreeClassifier(DecisionTree, ClassifierMixin):
    """Classify a multi-labeled dataset with a decision tree.
    
    Example
    -------
    
    >>> from scikits.learn.datasets import load_iris
    >>> from scikits.learn.cross_val import StratifiedKFold
    >>> data = load_iris()
    >>> train, test = StratifiedKFold(data.target, 10).next()
    >>> tree = DecisionTreeClassifier()
    >>> tree.fit(data.data[train], data.target[train])
    >>> tree.predict(data.data[test])

    """
    def __init__(self, criterion=information_gain, min_split=4, subsample=None, R=None):
        super(DecisionTreeClassifier, self).__init__(criterion=criterion)
        
    def fit(self, X, y, normalisedlabels=False, weights=None):
        if not normalisedlabels:
            labels, self.classes = normaliselabels(y)
        return super(DecisionTreeClassifier, self).fit(X, labels, weights)
    
    def predict(self, X):
        values = super(DecisionTreeClassifier, self).predict(X)
        return np.round(values)

class DecisionTreeRegressor(DecisionTree, RegressorMixin):
    raise NotImplementedError('Coming soon')

