# -*- coding: utf-8 -*-
# Copyright (C) 2008-2011, Luis Pedro Coelho <luis@luispedro.org>
# vim: set ts=4 sts=4 sw=4 expandtab smartindent:
#
# License: MIT. See COPYING.MIT file in the milk distribution

from __future__ import division
from collections import defaultdict
import numpy as np
from .classifier import normaliselabels

__all__ = [
    'zscore',
    'zscore_normalise',
    'interval_normalise',
    'chkfinite',
    'sample_to_2min',
    'normaliselabels'
]

def zscore(features):
    """
    features = zscore(features)

    Returns a copy of features which has been normalised to zscores 
    """
    mu = features.mean(0)
    sigma = np.std(features,0)
    sigma[sigma == 0] = 1
    return (features - mu) / sigma

class subtract_divide_model(object):
    def __init__(self, shift, factor):
        factor[factor == 0] = 1 # This makes the division a null op.

        self.shift = shift
        self.factor = factor

    def apply(self, features):
        return (features - self.shift)/self.factor

    def __repr__(self):
        return 'subtract_divide_model(%s, %s)' % (self.shift, self.factor)

class zscore_normalise(object):
    '''
    Normalise to z-scores

    A preprocessor that normalises features to z scores.
    '''

    def train(self, features, labels, normalisedlabels=False):
        shift = features.mean(0)
        factor = np.std(features,0)
        return subtract_divide_model(shift, factor)

class interval_normalise(object):
    '''
    Linearly scale to the interval [-1,1] (per libsvm recommendation)

    '''
    def train(self, features, labels, normalisedlabels=False):
        ptp = features.ptp(0)
        shift = features.min(0) + ptp/2.
        factor = ptp/2.
        return subtract_divide_model(shift, factor)

    def __repr__(self):
        return 'interval_normalise()'


def sample_to_2min(labels):
    '''
    selected = sample_to_2min(labels)

    Select examples so that the ratio of size of the largest
    class to the smallest class is at most two (i.e.,
        min_label_count = min { (labels == L).sum() | for L in set(labels) }
        for L' in set(labels):
            assert (labels == L').sum() <= 2 * min_label_count
    )

    Parameters
    ----------
        * labels: sequence of labels

    Output
    ------
        * selected: a Boolean numpy.ndarray
    '''
    counts = defaultdict(int)
    for n in labels:
        counts[n] += 1

    labels = np.asanyarray(labels)
    max_entries = np.min(counts.values())*2
    selected = np.zeros(len(labels), bool)
    for c in counts.iterkeys():
        p, = np.where(labels == c)
        p = p[:max_entries]
        selected[p] = 1
    return selected



class chkfinite(object):
    '''
    Fill NaN & Inf values

    Replaces NaN & Inf values with zeros.
    '''
    def __init__(self):
        pass

    def train(self, features, labels, normalisedlabels=False):
        return self

    def apply(self, features):
        nans = np.isnan(features) + np.isinf(features)
        if nans.any():
            features = features.copy()
            features[nans] = 0
        return features

    def __repr__(self):
        return 'chkfinite()'

