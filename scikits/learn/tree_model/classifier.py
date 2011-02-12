# -*- coding: utf-8 -*-
# Copyright (C) 2008-2011, Luis Pedro Coelho <luis@luispedro.org>
# vim: set ts=4 sts=4 sw=4 expandtab smartindent:
#
# Permission is hereby granted, free of charge, to any person obtaining a copy
#  of this software and associated documentation files (the "Software"), to deal
#  in the Software without restriction, including without limitation the rights
#  to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
#  copies of the Software, and to permit persons to whom the Software is
#  furnished to do so, subject to the following conditions:
#
# The above copyright notice and this permission notice shall be included in
#  all copies or substantial portions of the Software.
#
# THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
#  IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
#  FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
#  AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
#  LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
#  OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN
#  THE SOFTWARE.

from __future__ import division
import numpy as np

__all__ = ['normaliselabels', 'ctransforms']

def normaliselabels(labels):
    '''
    normalised, names = normaliselabels(labels)

    Normalises the labels to be integers from 0 through N-1

    `normalised` is a np.array, while `names` is a list mapping the indices to
    the old names.

    Parameters
    ----------
    labels : any iterable of labels

    Returns
    ------
    normalised : a numpy ndarray of integers 0 .. N-1
    names : list of label names
    '''
    names = sorted(set(labels))
    normalised = map(names.index, labels)
    return np.array(normalised), names

class threshold_model(object):
    '''
    threshold_model

    Attributes
    ----------
    threshold : float
        threshold value
    '''
    def __init__(self, threshold=.5):
        self.threshold = .5

    def apply(self, f):
        return f >= self.threshold

class fixed_threshold_learner(object):
    def __init__(self, threshold=.5):
        self.threshold = threshold
    def train(self, features, labels, **kwargs):
        return threshold_model(self.threshold)


class ctransforms_model(object):
    '''
    model = ctransforms_model(models)

    A model that consists of a series of transformations.

    See Also
    --------
      ctransforms
    '''
    def __init__(self, models):
        self.models = models
    def apply(self,features):
        for T in self.models:
            features = T.apply(features)
        return features

class ctransforms(object):
    '''
    ctransf = ctransforms(c0, c1, c2, ...)

    Concatenate transforms.
    '''
    def __init__(self,*args):
        self.transforms = args


    def train(self, features, labels, normalisedlabels=False):
        models = []
        model = None
        for T in self.transforms:
            if model is not None:
                features = np.array([model.apply(f) for f in features])
            model = T.train(features, labels, normalisedlabels=normalisedlabels)
            models.append(model)
        return ctransforms_model(models)


    def set_option(self, opt, val):
        idx, opt = opt
        self.transforms[idx].set_option(opt,val)

