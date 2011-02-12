# -*- coding: utf-8 -*-
# Copyright (C) 2010-2011, Luis Pedro Coelho <lpc@cmu.edu>
# vim: set ts=4 sts=4 sw=4 expandtab smartindent:
#
# License: MIT. See COPYING.MIT file in the milk distribution

'''
Random Forest
-------------

Main elements
-------------

rf_learner : A learner object
'''

from __future__ import division
import numpy as np
import milk.supervised.tree
from .normalise import normaliselabels

__all__ = [
    'rf_learner',
    ]

def _sample(features, labels, n, R):
    '''
    features', labels' = _sample(features, labels, n, R)

    Sample n element from (features,labels)

    Parameters
    ----------
    features : sequence
    labels : sequence
        Same size as labels
    n : integer
    R : random object

    Returns
    -------
    features' : sequence
    labels' : sequence
    '''

    N = len(features)
    sfeatures = []
    slabels = []
    for i in xrange(n):
        idx = R.randint(N)
        sfeatures.append(features[idx])
        slabels.append(labels[idx])
    return np.array(sfeatures), np.array(slabels)

class rf_model(object):
    def __init__(self, forest, names):
        self.forest = forest
        self.names = names

    def apply(self, features):
        rf = len(self.forest)
        votes = sum(t.apply(features) for t in self.forest)
        return (votes > (rf//2))
        

class rf_learner(object):
    '''
    Random Forest Learner

    learner = rf_learner(rf=101, frac=.7)

    Attributes
    ----------
    rf : integer, optional
        Nr of trees to learn (default: 101)
    frac : float, optional
        Sample fraction
    '''
    def __init__(self, rf=101, frac=.7):
        self.rf = rf
        self.frac = frac

    def train(self, features, labels, normalisedlabels=False, names=None, **kwargs):
        N,M = features.shape
        m = int(self.frac*M)
        n = int(self.frac*M)
        R = np.random
        tree = milk.supervised.tree.tree_learner()
        forest = []
        if not normalisedlabels:
            labels,names = normaliselabels(labels)
        elif names is None:
            names = (0,1)
        for i in xrange(self.rf):
            forest.append(
                    tree.train(*_sample(features, labels, n, R),
                               **{'normalisedlabels' : True})) # This syntax is necessary for Python 2.5
        return rf_model(forest, names)


