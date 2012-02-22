"""
Photometric redshifts using k-nearest neighbors

Usage:

% python photoz_nearest_neighbors.py data/sdss_photoz/

First argument is the data directory where the data has been
downloaded.  Note that you must first run the 'fetch_data.py' script
in data/sdss_phtoz/
"""
import sys, os

import numpy as np
import pylab as pl

from sklearn.neighbors import KNeighborsRegressor

try:
    data_dir = sys.argv[1]
    data_file = os.path.join(data_dir, 'sdss_galaxy_colors.npy')
    data = np.load(data_file)
except:
    print __doc__
    exit()

n_neighbors = 1

N = len(data)

# shuffle data
np.random.seed(0)
np.random.shuffle(data)

# put colors in a matrix
X = np.zeros((N, 4))
X[:, 0] = data['u'] - data['g']
X[:, 1] = data['g'] - data['r']
X[:, 2] = data['r'] - data['i']
X[:, 3] = data['i'] - data['z']
z = data['redshift']

# divide into training and testing data
Ntrain = N/2
Xtrain = X[:Ntrain]
ztrain = z[:Ntrain]

Xtest = X[Ntrain:]
ztest = z[Ntrain:]

knn = KNeighborsRegressor(n_neighbors, weights='uniform')
zpred = knn.fit(Xtrain, ztrain).predict(Xtest)

axis_lim = np.array([-0.1, 2.5])

rms = np.sqrt(np.mean((ztest - zpred) ** 2))
print rms
print len(ztest)
print np.sum(abs(ztest - zpred) > 1)

ax = pl.axes()
pl.scatter(ztest, zpred, c='k', lw=0, s=4)
pl.plot(axis_lim, axis_lim, '--k')
pl.plot(axis_lim, axis_lim + rms, ':k')
pl.plot(axis_lim, axis_lim - rms, ':k')
pl.xlim(axis_lim)
pl.ylim(axis_lim)

pl.text(0.99, 0.02, "RMS error = %.2g" % rms,
        ha='right', va='bottom', transform=ax.transAxes,
        bbox=dict(ec='w', fc='w'), fontsize=16)

pl.title('Photo-z: Nearest Neigbor Regression')
pl.xlabel(r'$\mathrm{z_{spec}}$', fontsize=14)
pl.ylabel(r'$\mathrm{z_{phot}}$', fontsize=14)
pl.show()



