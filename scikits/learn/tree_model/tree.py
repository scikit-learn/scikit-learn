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
from .classifier import normaliselabels

__all__ = [
    'tree_learner',
    'stump_learner',
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


from ._tree import set_entropy
from ._tree import information_gain as _information_gain
def information_gain(labels0, labels1, include_entropy=False):
    '''
    ig = information_gain(labels0, labels1, include_entropy=False)

    Information Gain
    See http://en.wikipedia.org/wiki/Information_gain_in_decision_trees

    The function calculated here does not include the original entropy unless
    you explicitly ask for it (by passing include_entropy=True)
    '''
    if include_entropy:
        return set_entropy(np.concatenate( (labels0, labels1) )) + \
                _information_gain(labels0, labels1)
    return _information_gain(labels0, labels1)


def z1_loss(labels0, labels1, weights0=None, weights1=None):
    '''
    z = z1_loss(labels0, labels1)
    z = z1_loss(labels0, labels1, weights0, weights1)

    zero-one loss split for tree learning
    '''
    def _acc(labels, weights):
        c = (labels.mean() > .5)
        if weights is not None:
            return -np.dot((labels != c), weights)
        return -np.sum(labels != c)
    return _acc(labels0, weights0) + _acc(labels1, weights1)



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
        if given, then, at each step, choose
    R : source of randomness, optional
        See `milk.util.get_pyrandom`
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
            raise ValueError('milk.supervised.tree.build_tree: `subsample` must be > 0.\nDid you mean to use None to signal no subsample?')
        from ..utils import get_pyrandom
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


class tree_learner(object):
    '''
    tree = tree_learner()
    model = tree.train(features, labels)
    model2 = tree.train(features, labels, weights=weights)
    predicted = model.apply(testfeatures)

    A decision tree classifier (currently, implements the greedy ID3
    algorithm without any pruning).

    Attributes
    ----------
    criterion : function, optional
        criterion to use for tree construction,
        this should be a function that receives a set of labels
        (default: information_gain).

    min_split : integer, optional
        minimum size to split on (default: 4).
    '''
    def __init__(self, criterion=information_gain, min_split=4, return_label=True, subsample=None, R=None):
        self.criterion = criterion
        self.min_split = min_split
        self.return_label = return_label
        self.subsample = subsample
        self.R = R

    def train(self, features, labels, normalisedlabels=False, weights=None):
        if not normalisedlabels:
            labels,names = normaliselabels(labels)
        tree = build_tree(features, labels, self.criterion, self.min_split, self.subsample, self.R, weights)
        return tree_model(tree, self.return_label)

tree_classifier = tree_learner

class tree_model(object):
    '''
    tree model
    '''
    def __init__(self, tree, return_label):
        self.tree = tree
        self.return_label = return_label

    def apply(self,feats):
        value = apply_tree(self.tree, feats)
        if self.return_label:
            return value > .5
        return value

class stump_model(object):
    def __init__(self, idx, cutoff, names):
        self.names = names
        self.idx = idx
        self.cutoff = cutoff

    def apply(self, f):
        return self.names[f[self.idx] > self.cutoff]

    def __repr__(self):
        return '<stump(%s, %s)>' % (self.idx, self.cutoff)

class stump_learner(object):
    def __init__(self):
        pass

    def train(self, features, labels, normalisedlabels=False, weights=None, **kwargs):
        if not normalisedlabels:
            labels,names = normaliselabels(labels)
        idx,cutoff = _split(features, labels, weights, z1_loss, subsample=None, R=None)
        return stump_model(idx, cutoff, names)
