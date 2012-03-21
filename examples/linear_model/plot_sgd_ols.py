"""
===============================
Ordinary Least Squares with SGD
===============================

Simple Ordinary Least Squares example with stochastic
gradient descent, we draw the linear least
squares solution for a random set of points in the plane.
"""
print __doc__

import pylab as pl

from sklearn.linear_model import SGDRegressor
from sklearn.datasets.samples_generator import make_regression

# this is our test set, it's just a straight line with some
# gaussian noise
X, Y = make_regression(n_samples=100, n_features=1, n_informative=1,\
                        random_state=0, noise=35)

# run the classifier
clf = SGDRegressor(alpha=0.1, n_iter=20)
clf.fit(X, Y)

# and plot the result
pl.scatter(X, Y, color='black')
pl.plot(X, clf.predict(X), color='blue', linewidth=3)
pl.show()
