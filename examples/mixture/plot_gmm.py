"""
=================================
Gaussian Mixture Model Ellipsoids
=================================

Plot the confidence ellipsoids of a mixture of two gaussians with EM
and variational dirichlet process.

Both models have access to five components with which to fit the
data. Note that the EM model will necessarily use all five components
while the DP model will effectively only use as many as are needed for
a good fit. This is a property of the Dirichlet Process prior.

This example doesn't show it, as we're in a low-dimensional space, but
another advantage of the dirichlet process model is that it can fit
full covariance matrices effectively even when there are less examples
per cluster than there are dimensions in the data, due to
regularization properties of the inference algorithm.
"""

import numpy as np
from scikits.learn import mixture
import itertools

import pylab as pl
import matplotlib as mpl

n, m = 300, 2

# generate random sample, two components
np.random.seed(0)
C = np.array([[0., -0.7], [3.5, .7]])
X = np.r_[np.dot(np.random.randn(n, 2), C),
          np.random.randn(n, 2) + np.array([3, 3])]


# fit a mixture of gaussians with EM using five components
clf = mixture.GMM(n_states=5, cvtype='diag')
clf.fit(X)

# fit a dirichlet process mixture of gaussians using five components
dpclf = mixture.DPGMM(n_states=5, cvtype='diag')
dpclf.fit(X)

color_iter = itertools.cycle (['r', 'g', 'b', 'c', 'm'])


for i,c in enumerate([clf, dpclf]):
    splot = pl.subplot(211+i, aspect='equal')
    Y_ = c.predict(X)
    for i, (mean, covar, color) in enumerate(zip(c.means, c.covars, color_iter)):
        v, w = np.linalg.eigh(covar)
        u = w[0] / np.linalg.norm(w[0])
        # as the DP will not use every component it has access to
        # unless it needs it, we shouldn't plot the redundant
        # components.
        if not sum(Y_ == i) > 1:
            continue
        pl.scatter(X[Y_==i, 0], X[Y_==i, 1], .8, color=color)
        angle = np.arctan(u[1]/u[0])
        angle = 180 * angle / np.pi # convert to degrees
        ell = mpl.patches.Ellipse (mean, v[0], v[1], 180 + angle, color=color)
        ell.set_clip_box(splot.bbox)
        ell.set_alpha(0.5)
        splot.add_artist(ell)

# Note that the GMM will use all components it has access to, while
# the dirichlet process model will only use as many are needed to
# explain the data
pl.show()

