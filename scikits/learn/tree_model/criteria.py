# -*- coding: utf-8 -*-
# Copyright (C) 2008-2011, Luis Pedro Coelho <luis@luispedro.org>
# vim: set ts=4 sts=4 sw=4 expandtab smartindent:
#
# License: MIT. See COPYING.MIT file in the milk distribution
'''
==============
Tree Splitters
==============

Criteria options for tree splitting

'''

from __future__ import division
import numpy as np
from .base import normaliselabels

__all__ = [
    'information_gain',
    'z1_loss'
    ]

from ._tree import set_entropy, gini_index

def information_gain(left_labels, right_labels, include_entropy=True):
    """
    http://en.wikipedia.org/wiki/Information_gain_in_decision_trees
    """
    clen = max(left_labels.max(), right_labels.max()) + 1
    N0 = np.prod(left_labels.shape)
    N1 = np.prod(right_labels.shape)
    N = N0 + N1
    H = - N0/N * set_entropy(left_labels, N0,  clen) - N1/N * set_entropy(right_labels, N1, clen)
    if include_entropy:
        H += set_entropy(np.concatenate( (left_labels, right_labels) ), N, clen)
    return H

def information_gain_gini(left_labels, right_labels, include_entropy=True):
    """
    http://en.wikipedia.org/wiki/Information_gain_in_decision_trees
    """
    clen = max(left_labels.max(), right_labels.max()) + 1
    N0 = np.prod(left_labels.shape)
    N1 = np.prod(right_labels.shape)
    N = N0 + N1
    H = - N0/N * gini_index(left_labels, N0,  clen) - N1/N * gini_index(right_labels, N1, clen)
    if include_entropy:
        H += gini_index(np.concatenate( (left_labels, right_labels) ), N, clen)
    return H

def information_gain_ratio(left_labels, right_labels):
    """
    http://en.wikipedia.org/wiki/Information_gain_in_decision_trees
    """
    IG = information_gain(left_labels, right_labels)
    N = np.prod(left_labels.shape) + np.prod(right_labels.shape)
    IV = -N*(1/N)*np.log(1/N) 
    return IG/IV

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
