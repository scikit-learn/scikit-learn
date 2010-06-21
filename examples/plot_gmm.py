"""
Simple Gaussian Mixture model plotting example

TODO: use the faithful dataset
"""

import numpy as np
from scikits.learn import gmm, shortcuts
from datetime import datetime
import itertools

import pylab as pl
import matplotlib as mpl

n, m = 442, 2

np.random.seed(0)

C = np.array([[0., -0.7], [3.5, .7]]) # rotation and scale matrix

X = np.r_[np.dot(np.random.randn(n, 2), C),
          np.random.randn(n, 2) + np.array([10,10])]

clf = gmm.GMM(2, cvtype='full')
clf.fit(X)

splot = pl.subplot(111, aspect='equal')

n_comp = len(clf.means) # number of components

color_iter = itertools.cycle (['r', 'g', 'b', 'c'])

Y_ = clf.predict(X)

for i, mean, covar, color in zip(range(n_comp), clf.means,
                                 clf.covars, color_iter):
    v, w = np.linalg.eigh(covar)
    u = w[0] / np.linalg.norm(w[0])
    pl.scatter(X[Y_==i, 0], X[Y_==i, 1], .8, color=color)
    angle = np.arctan(u[1]/u[0])
    angle = 180 * angle / np.pi # convert to degrees
    ell = mpl.patches.Ellipse (mean, v[0], v[1], 180 + angle, color=color) #, angle)
    ell.set_clip_box(splot.bbox)
    ell.set_alpha(0.5)
    splot.add_artist(ell)

pl.show()

