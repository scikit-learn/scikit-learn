"""
This package collects helper functions and classes that span multiple
modules of scikits.learn. In other words, these functions introduce
controlled coupling for convenience's sake.
"""

import numpy as np
import itertools

def plot_gmixture (X, means, covars, r=1.):
    """
    Plot a mixture of gaussians.

    requires matplotlib

    Work in progress
    
    """
    import matplotlib as mpl
    import pylab as pl

    assert X.shape[1] == 2
    splot = pl.subplot(111, aspect='equal')

    pl.scatter(X[:, 0], X[:, 1], .8)

    color_iter = itertools.cycle (['r', 'g', 'b', 'c'])

    for mean, covar, color in zip(means, covars, color_iter):
        v, w = np.linalg.eigh(covar) # todo: not general!
        v *= r
        angle = np.arccos(w[0, 0]/w[0, 1])
        ell = mpl.patches.Ellipse (mean, v[0], v[1], color=color) #, angle)
        ell.set_clip_box(splot.bbox)
        ell.set_alpha(0.5)
        splot.add_artist(ell)


    pl.show()

    

    
