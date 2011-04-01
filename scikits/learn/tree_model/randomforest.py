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
from .base import normaliselabels
from .tree import DecisionTree

__all__ = [
    'RandomForest',
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

class RandomForestBase(object):
    '''Random Forest base class

    Attributes
    ----------
    
    criterion_func : function handle
        criterion_function to use
    num_trees : integer, optional
        Nr of trees to learn (default: 101)
    sample_fraction : float, optional
        Sample fraction
        
    '''
    
    def __init__(self, num_trees=101, sample_fraction=.7):
        self.num_trees = num_trees
        self.sample_fraction = sample_fraction
        self.forest = None
        self.names = None

    def fit(self, features, labels, normalisedlabels=False, names=None, **kwargs):
        N,M = features.shape
        m = int(self.sample_fraction*M)
        n = int(self.sample_fraction*M)
        R = np.random
        forest = []
        if not normalisedlabels:
            labels,names = normaliselabels(labels)
        elif names is None:
            names = (0,1)
        for i in xrange(self.num_trees):
            tree = DecisionTree()
            tree.fit(*_sample(features, labels, n, R),
                      **{'normalisedlabels' : True}) # This syntax is necessary for Python 2.5
            forest.append(tree)
        self.forest = forest
        self.names = names
        
    def predict(self, features):
        num_trees = len(self.forest)
        votes = sum(t.predict(features) for t in self.forest)
        return (votes > (num_trees//2))

class RandomForest(RandomForestBase):
    pass
