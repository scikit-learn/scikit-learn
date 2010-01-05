#! /usr/bin/env python
# Last Change: Mon Jul 09 05:00 PM 2007 J

# Various utilities for examples 

import numpy as N
from numpy.testing import set_package_path, restore_path

from scikits.learn.datasets import oldfaithful, pendigits, iris

def get_faithful():
    """Return faithful data as a nx2 array, first column being duration, second
    being waiting time."""
    # Load faithful data, convert waiting into integer, remove L, M and S data
    data = oldfaithful.load()
    tmp1 = []
    tmp2 = []
    for i in data:
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
    data = pendigits.load()
    return data['training']['x'], data['training']['y']

def scale(data, mode = 'sym'):
    """Linearly scale data in place such as each col is in the range [0..1].

    Returns the translation factor t and scaling factor s. You can retrieve
    the original values with data = s * scaled + t."""
    n = N.min(data, 0)
    m = N.max(data, 0)
    if mode == 'sym':
        t = n + 0.5 * (m - n)
        s = 0.5 * (m - n)
    elif mode == 'right':
        t = n
        s = m - n
    else:
        raise ValueError("Mode %s not recognized" % mode)
    
    data -= t
    data /= s
    return t, s

if __name__ == '__main__':
    a = N.random.randn(10, 2)
    b = a.copy()
    scale(a)
    print a
    scale(b, 'right')
    print b
