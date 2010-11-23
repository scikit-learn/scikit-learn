# -*- coding: utf-8 -*-

import numpy as np

def dist2hd(x,y):
    """
    Generates a distance matrix

    Parameters
    ----------
    x : array_like

    y : array_like

    Returns
    -------
    The distance matrix between x and y
    """
    d = np.zeros((x.shape[0],y.shape[0]),dtype=x.dtype)
    for i in xrange(x.shape[1]):
        diff2 = x[:,i,None] - y[:,i]
        diff2 **= 2
        d += diff2
    np.sqrt(d,d)
    return d
 
