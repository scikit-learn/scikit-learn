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
from .tree import DecisionTreeClassifier
from ..base import BaseEstimator, ClassifierMixin, RegressorMixin

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

    N, _ = features.shape
    sfeatures = []
    slabels = []
    for i in xrange(n):
        idx = R.randint(N)
        sfeatures.append(features[idx])
        slabels.append(labels[idx])
    return np.array(sfeatures), np.array(slabels)

class RandomForestBase(BaseEstimator):
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
    
    def __init__(self, num_trees=101, sample_fraction=.5):
        self.num_trees = num_trees
        self.sample_fraction = sample_fraction
        self.forest = None

    def fit(self, X, y, normalisedlabels=False):
        observations, features = X.shape
        n = int(self.sample_fraction*observations)
        R = np.random
        forest = []
        labels = y
        self.classes = np.unique(y)
        if not normalisedlabels:
            labels, names = normaliselabels(y)
            self.classes = names
        for i in xrange(self.num_trees):
            tree = DecisionTreeClassifier()
            tree.fit(*_sample(X, labels, n, R),
                      **{'normalisedlabels' : True}) # This syntax is necessary for Python 2.5
            forest.append(tree)
        self.forest = forest
        return self
        
    def predict(self, X):
        n_observations,_ = X.shape
        values = np.zeros(n_observations)
        for idx, f in enumerate(X):
            out = [t.predict(f[None,:]).item() for t in self.forest]
            votes = np.bincount(out)
            values[idx] = votes.argmax()
        return values

class RandomForest(RandomForestBase, ClassifierMixin):
    pass
