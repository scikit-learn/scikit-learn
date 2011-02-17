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

from ._tree import set_entropy
def information_gain(labels0, labels1, include_entropy=False):
    clen = max(labels0.max(), labels1.max()) + 1
    N0 = np.prod(labels0.shape)
    N1 = np.prod(labels1.shape)
    N = N0 + N1
    H = - N0/N * set_entropy(labels0, N0,  clen) - N1/N * set_entropy(labels1, N1, clen)
    if include_entropy:
        H += set_entropy(np.concatenate( (labels0, labels1) ), N, clen)
    return H

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
