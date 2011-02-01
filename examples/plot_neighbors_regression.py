"""
==============================
k-Nearest Neighbors regression
==============================

Demonstrate the resolution of a regression problem
using a k-Nearest Neighbor and the interpolation of the
target using barycenter computation.

"""
print __doc__

###############################################################################
# Generate sample data
import numpy as np

np.random.seed(0)
X = np.sort(5*np.random.rand(40, 1), axis=0)
T = np.linspace(0, 5, 500)[:, np.newaxis]
y = np.sin(X).ravel()

# Add noise to targets
y[::5] += 1*(0.5 - np.random.rand(8))

###############################################################################
# Fit regression model

from scikits.learn import neighbors

knn_barycenter = neighbors.NeighborsBarycenter(n_neighbors=5)
y_ = knn_barycenter.fit(X, y).predict(T)

###############################################################################
# look at the results
import pylab as pl
pl.scatter(X, y, c='k', label='data')
pl.hold('on')
pl.plot(T, y_, c='g', label='k-NN prediction')
pl.xlabel('data')
pl.ylabel('target')
pl.legend()
pl.title('k-NN Regression')
pl.show()

