"""
This package collects helper functions and classes that span multiple
modules of scikits.learn. In other words, these functions introduce
controlled coupling for convenience's sake.
"""

import numpy as np
import itertools

def plot_gmixture (X, means, covars, clf=None, r=1.):
    """
    Plot a mixture of gaussians.

    requires matplotlib

    Parameters
    ----------

    X : array, [nsamples, nfeatures]
        The sample to plot

    mean : array,
        vector of means.

    covars : array
        vector of covariances

    clf : classifier, optional
        If given, will plot also the sample in different colors for
        depending on their classification results

    r : radious of plotted ellipsis
        If not given, r = 1

    Examples
    --------
    See examples/plot_gmm.py for a complete example.
    """
    import matplotlib as mpl
    import pylab as pl

    assert X.shape[1] == 2
    splot = pl.subplot(111, aspect='equal')

    n_comp = len(means) # number of components

    color_iter = itertools.cycle (['r', 'g', 'b', 'c'])

    if clf is None:
        pl.scatter(X[:, 0], X[:, 1], .8)
    else:
        Y_ = clf.predict(X)

    for i, mean, covar, color in zip(range(n_comp), means,
                                     covars, color_iter):
        v, w = np.linalg.eigh(covar)
        v *= r
        u = w[0] / np.linalg.norm(w[0])
        if clf is not None:
            pl.scatter(X[Y_==i, 0], X[Y_==i, 1], .8, color=color)
        angle = np.arctan(u[1]/u[0])
        angle = 180 * angle / np.pi # convert to degrees
        ell = mpl.patches.Ellipse (mean, v[0], v[1], 180 + angle, color=color) #, angle)
        ell.set_clip_box(splot.bbox)
        ell.set_alpha(0.5)
        splot.add_artist(ell)


    pl.show()

    

    
