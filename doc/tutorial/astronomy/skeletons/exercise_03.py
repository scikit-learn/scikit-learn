"""
Astronomy Tutorial: exercise 3

Dimensionality reduction of stellar spectra

Usage: python exercise_03.py datadir [-m method] [-k n_neigbors]
                                     [-n norm_type] [-N n_samples]
                                     [-s]

  - datadir is $TUTORIAL_DIR/data/sdss_photoz
    This directory should contain the file sdss_photoz.npy

  - method is one of [pca | lle | mlle | isomap].  If not specified,
    PCA will be performed
    
  - n_neighbors is an integer number of neighbors to use with manifold methods

  - norm_type is one of [none | l1 | l2].  It specifies how the data should
    be normalized.

  - n_samples is the number of samples used for the projection.  Only 1000
    of the 4000 samples are used by default.

  - specifying -s shuffles the data.  This can help test for stability of
    the reconstruction.

Description:
In this tutorial, we explore manifold learning techniques to visualize 4000
SDSS spectral data.  This is a much more exploratory exercise than the previous
two.  The goal is to determine how to best visualize this high-dimensional
space.  You will implement PCA, LLE, Modified LLE, and Isomap, for various
data normalizations.  The goal is to find the best visualization of the
data, where "best" in this case is a qualitative measure of how well the
different classes of points are separated in the projected space.

To make experimentation more streamlined

There are several places in this file with code to be filled-in as part of
the exercise.  Each of these is labeled TODO below.
"""

import os, sys
import numpy as np

import pylab as pl
from matplotlib import ticker

from sklearn import preprocessing
from sklearn.decomposition import RandomizedPCA
from sklearn.manifold import LocallyLinearEmbedding, Isomap

#----------------------------------------------------------------------
# set up command-line option parser
from optparse import OptionParser
parser = OptionParser(usage=__doc__,
                      version="%prog 1.0")
parser.add_option("-m", "--method",
                  dest="method",
                  default='pca',
                  help="Specify method to use: [pca | lle | mlle | isomap]")

parser.add_option("-k", "--neighbors",
                  dest="n_neighbors",
                  type="int",
                  default=15,
                  help='Specify number of neighbors for manifold learning')

parser.add_option("-N", "--normalization",
                  dest="normalization",
                  default="none",
                  help="Specify normalization: [none | l1 | l2]")

parser.add_option("-n", "--n_samples",
                  dest="n_samples",
                  type="int",
                  default=1000,
                  help="Specify number of samples to use, up to 4000 (default 1000)")

parser.add_option("-s", "--shuffle",
                  dest="shuffle",
                  action="store_true",
                  default=False,
                  help="shuffle the data")


options, args = parser.parse_args()

if len(args) == 0:
    parser.error("Must specify a data directory")
elif len(args) > 1:
    parser.error("Must specify a single data directory")

datadir = args[0]

print "data directory: %s" % datadir
print " method = %s" % options.method
print " n_neighbors = %i" % options.n_neighbors
print " normalization = %s" % options.normalization
print " n_samples: %i" % options.n_samples
print " shuffle: %s" % options.shuffle


def three_component_plot(c1, c2, c3, color, labels):
    pl.figure(figsize=(8,8))
    kwargs = dict(s=4, lw=0, c=color, vmin=2, vmax=6)
    ax1 = pl.subplot(221)
    pl.scatter(c1, c2, **kwargs)
    pl.ylabel('component 2')

    ax2 = pl.subplot(223, sharex=ax1)
    pl.scatter(c1, c3, **kwargs)
    pl.xlabel('component 1')
    pl.ylabel('component 3')

    ax3 = pl.subplot(224, sharey=ax2)
    pl.scatter(c2, c3, **kwargs)
    pl.xlabel('component 2')

    for ax in (ax1, ax2, ax3):
        ax.xaxis.set_major_formatter(ticker.NullFormatter())
        ax.yaxis.set_major_formatter(ticker.NullFormatter())

    pl.subplots_adjust(hspace=0.05, wspace=0.05)

    format = ticker.FuncFormatter(lambda i, *args: labels[i])
    pl.colorbar(ticks = range(2, 7), format=format,
                cax = pl.axes((0.52, 0.51, 0.02, 0.39)))
    pl.clim(1.5, 6.5)


#----------------------------------------------------------------------
# Load data files
data = np.load(os.path.join(datadir, 'spec4000_corrected.npz'))

X = data['X']
y = data['y']
labels = data['labels']

if options.shuffle:
    i = np.arange(y.shape[0], dtype=int)
    np.random.shuffle(i)
    X = X[i]
    y = y[i]

#----------------------------------------------------------------------
# truncate the data for experimentation
#
#  There are 4000 points, which can take a long time to run.  By default,
#  it is truncated to 1000 samples.  This can be changed using the -n
#  command-line argument.

X = X[:options.n_samples]
y = y[:options.n_samples]

#----------------------------------------------------------------------
# Normalization: 
#  
#  The results of the dimensionality reduction can depend heavily on the
#  data normalization.  These can be commented or un-commented to try
#  l1 or l2 normalization.

if options.normalization.lower() == 'none':
    pass
elif options.normalization.lower() == 'l2':
    X = preprocessing.normalize(X, 'l2')
elif options.normalization.lower() == 'l1':
    X = preprocessing.normalize(X, 'l1')
else:
    raise ValueError("Unrecognized normalization: '%s'" % options.normalization)

#======================================================================
# TODO: compute X_proj for each method.
#   In each of the below cases, you should compute a projection of the
#   data and store that projection in the matrix X_proj.
#   X_proj should have the same number of rows as X, and should have
#   at least 3 features.

X_proj = None

if options.method == 'pca':
    print "Performing PCA"
    # TODO:  compute a RandomizedPCA projection of X with n_components >= 3

elif options.method == 'lle':
    print "Performing LLE"
    # TODO:  compute LLE on X with method='standard', and out_dim >= 3


elif options.method == 'mlle':
    print "Performing MLLE"
    # TODO:  compute LLE on X with method='modified' and out_dim >= 3

elif options.method == 'isomap':
    print "Performing Isomap"
    # TODO:  compute Isomap on X with out_dim >= 3

else:
    raise ValueError("Unrecognized method: '%s'" % options.method)

three_component_plot(X_proj[:, 0], X_proj[:, 1], X_proj[:, 2], y, labels)
pl.show()

