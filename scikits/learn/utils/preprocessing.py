#! /usr/bin/env python
# Last Change: Mon Jul 09 08:00 PM 2007 J

# Various utilities for examples 

import numpy as N

"""Different tools for pre processing, like whitening or scaling data."""

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

def whiten():
    """Whiten data."""
    raise NotImplementedError("whitening not implemented yet")
