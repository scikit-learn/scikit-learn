"""
This package collects helper functions and classes that span multiple
levels of scikits.learn. In other words, these functions introduce
controlled coupling for convenience's sake.
"""

import numpy as np

def plot_gmixture (X, sigma, clf):
    """
    requires matplotlib

    Work in progress
    
    k : number of components
    
    TODO: better name. Can this be generalized to arbitrary mixture models ?
    """
    import matplotlib
    import pylab as pl

    assert X.shape[1] == 2
    a = pl.subplot(111, aspect='equal')

    for mean, covar in zip(clf.means, clf.covars):
        v, w = np.linalg.eigh(np.diag(covar)) # todo: not general!
        angle = np.arccos(w[0, 0]/w[0, 1])
        ell = matplotlib.patches.Ellipse (mean, v[0], v[1]) #, angle)
        ell.set_clip_box(a.bbox)
        ell.set_alpha(0.5)
        a.add_artist(ell)


    pl.scatter(X[:, 0], X[:, 1], .8)
    pl.show()

    

    
