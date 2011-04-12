"""
=================================
DP Mixture Model Ellipsoids
=================================

Plot the covariance ellipsoids of a dirichlet process mixture of two
gaussians for varying values of the alpha parameter.

Note that we generate the data from two components, which are
correctly recovered by the dirichlet process, even though the
approximating distribution is truncated at five components.
"""

import numpy as np
from scikits.learn import mixture
import itertools

import pylab as pl
import matplotlib as mpl

n, m = 200, 2

# generate random sample, two components
np.random.seed(0)
C = np.array([[0., -0.7], [3.5, .7]])
X = np.r_[np.dot(np.random.randn(n, 2), C),
          np.random.randn(n, 2) + np.array([3, 3])]

for p, alpha in enumerate([0.01, 1.]):
    # fit a give-component dirichlet process mixture model on the
    # data.
    clf = mixture.DPGMM(n_states=5, cvtype='diag', alpha=alpha)
    clf.fit(X)

    splot = pl.subplot(311 + p, aspect='equal')
    color_iter = itertools.cycle(['r', 'g', 'b', 'c', 'm', 'y'])

    Y_ = clf.predict(X)

    for i, (mean, covar, color) in enumerate(zip(clf.means,
                                                 clf.covars,
                                                 color_iter)):
        v, w = np.linalg.eigh(covar)
        u = w[0] / np.linalg.norm(w[0])
        # as the DP will not use every component it has access to
        # unless it needs it, we shouldn't plot the redundant
        # components.
        if not sum(Y_ == i) > 1:
            continue
        pl.scatter(X[Y_ == i, 0], X[Y_ == i, 1], .8, color=color)
        angle = np.arctan(u[1] / u[0])
        angle = 180 * angle / np.pi  # convert to degrees
        ell = mpl.patches.Ellipse(mean, v[0], v[1], 180 + angle, color=color)
        ell.set_clip_box(splot.bbox)
        ell.set_alpha(0.5)
        splot.add_artist(ell)

pl.show()
