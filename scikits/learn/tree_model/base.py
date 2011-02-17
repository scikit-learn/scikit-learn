# -*- coding: utf-8 -*-
# Copyright (C) 2008-2010, Luis Pedro Coelho <lpc@cmu.edu>
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
import random

__all__ = ['normaliselabels', 'get_nprandom', 'get_pyrandom']

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

def get_nprandom(R):
    '''
    R' = get_nprandom(R)

    Returns a numpy.RandomState from R

    Parameters
    ----------
    R : can be one of:
        None          : Returns the default numpy global state
        integer       : Uses it as a seed for constructing a new random generator
        RandomState   : returns R

    Returns
    -------
    R' : np.RandomState
    '''
    if R is None:
        return np.random.mtrand._rand
    if type(R) == int:
        return np.random.RandomState(R)
    if type(R) is random.Random:
        return np.random.RandomState(R.randint(0, 2**30))
    if type(R) is np.random.RandomState:
        return R
    raise TypeError,"get_nprandom() does not know how to handle type %s." % type(R)

def get_pyrandom(R):
    '''
    R = get_pyrandom(R)

    Returns a random.Random object based on R

    Parameters
    ----------
    R : can be one of:
        None          : Returns the default numpy global state
        integer       : Uses it as a seed for constructing a new random generator
        RandomState   : returns R

    Returns
    -------
    R' : random.Random
    '''
    if R is None:
        return random.seed.im_self
    if type(R) is int:
        return random.Random(R)
    if type(R) is np.random.RandomState:
        return random.Random(R.randint(2**30))
    if type(R) is random.Random:
        return R
    raise TypeError,"get_pyrandom() does not know how to handle type %s." % type(R)
