"""
===============================
Ordinary Least Squares with SGD
===============================

Simple Ordinary Least Squares example with stochastic
gradient descent, we draw the linear least
squares solution for a random set of points in the plane.
"""
print __doc__

import numpy as np
import pylab as pl

from sklearn.linear_model import SGDRegressor

# this is our test set, it's just a straight line with some
# gaussian noise
xmin, xmax = -5, 5
n_samples = 100
X = [[i] for i in np.linspace(xmin, xmax, n_samples)]
Y = 2 + 0.5 * np.linspace(xmin, xmax, n_samples) \
      + np.random.randn(n_samples, 1).ravel()

# run the classifier
clf = SGDRegressor(alpha=0.1, n_iter=20)
clf.fit(X, Y)

# and plot the result
pl.scatter(X, Y, color='black')
pl.plot(X, clf.predict(X), color='blue', linewidth=3)
pl.show()
