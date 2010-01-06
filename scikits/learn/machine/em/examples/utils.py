#! /usr/bin/env python
# Last Change: Sun Jul 22 03:00 PM 2007 J

# Various utilities for examples 

import numpy as N

from scikits.learn.datasets import oldfaithful, pendigits, iris

def get_faithful():
    """Return faithful data as a nx2 array, first column being duration, second
    being waiting time."""
    # Load faithful data, convert waiting into integer, remove L, M and S data
    data = oldfaithful.load()
    tmp1 = []
    tmp2 = []
    for i in data['data']:
        if not (i[0] == 'L' or i[0] == 'M' or i[0] == 'S'):
            tmp1.append(i[0])
            tmp2.append(i[1])
            
    waiting = N.array([int(i) for i in tmp1], dtype = N.float)
    duration = N.array([i for i in tmp2], dtype = N.float)

    waiting = waiting[:, N.newaxis]
    duration = duration[:, N.newaxis]

    return N.concatenate((waiting, duration), 1)

def get_pendigits():
    """Return faithful data as a nx2 array, first column being duration, second
    being waiting time."""
    # Load faithful data, convert waiting into integer, remove L, M and S data
    data = pendigits.training.load()
    return data['data']['x0'], data['data']['y0']

def scale(data):
    """ Scale data such as each col is in the range [0..1].

    Note: inplace."""
    n = N.min(data, 0)
    m = N.max(data, 0)

    data -= n
    data /= (m-n)
    return data

