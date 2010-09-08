"""
==========================
FastICA on 2D point clouds
==========================

XXX: Still buggy !!!

"""

import numpy as np
import pylab as pl
from scipy import stats

from scikits.learn.pca import PCA
from scikits.learn.fastica import FastICA

###############################################################################
# Generate sample data
S = np.c_[np.random.normal(scale=0.5, size=(2, 10000)),
          np.random.normal(scale=4, size=(2, 10000))]
np.random.shuffle(S)

def pdf(x):
    e = np.exp(1.)
    return 0.5*(stats.norm(scale=0.5/e).pdf(x) + stats.norm(scale=4/e).pdf(x))

# Standardize data
S /= S.std(axis=1)[:,np.newaxis]
density = pdf(S[0]) * pdf(S[1])

# Mix data
A = [[1, 1], [0, 2]] # Mixing matrix

X = np.dot(A, S) # Generate observations
X /= np.std(X)

pca = PCA()
S_pca_ = pca.fit(X.T).transform(X.T).T

ica = FastICA()
S_ica_ = ica.fit(X).transform(X) # Estimate the sources

S_ica_ /= S_ica_.std(axis=1)[:,np.newaxis]

###############################################################################
# Plot results

def plot_samples(S, density, axis_list=None):
    pl.scatter(S[0], S[1], s=1, marker='o', c=np.sqrt(density), linewidths=0,
                                            zorder=10)
    if axis_list is not None:
        colors = [(0, 0.6, 0), (0.6, 0, 0)]
        for color, axis in zip(colors, axis_list):
            axis /= axis.std()
            x_axis, y_axis = axis
            # Trick to get legend to work
            pl.plot(0.1*x_axis, 0.1*y_axis, linewidth=2, color=color)
            # pl.quiver(x_axis, y_axis, x_axis, y_axis, zorder=11, width=0.01,
            pl.quiver(0, 0, x_axis, y_axis, zorder=11, width=0.01,
                        scale=6, color=color)

    pl.hlines(0, -3, 3)
    pl.vlines(0, -3, 3)
    pl.xlim(-3, 3)
    pl.ylim(-3, 3)
    pl.xlabel('$x$')
    pl.ylabel('$y$')

pl.close('all')
pl.subplot(2, 2, 1)
plot_samples(S, density)
pl.title('True Independant Sources')

axis_list = [pca.components_, ica.get_mixing_matrix()]
# axis_list = [pca.components_, np.array(A)]
pl.subplot(2, 2, 2)
plot_samples(X, density, axis_list=axis_list)
pl.legend(['PCA', 'ICA'], loc='upper left')
pl.title('Observations')

pl.subplot(2, 2, 3)
plot_samples(S_pca_, density)
pl.title('PCA scores')

pl.subplot(2, 2, 4)
plot_samples(S_ica_, density)
pl.title('ICA sources')

pl.subplots_adjust(0.09, 0.04, 0.94, 0.94, 0.26, 0.26)

pl.show()