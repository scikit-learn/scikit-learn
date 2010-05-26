"""
Simple Ordinary Least Squares example.

We draw the linear least squares solution for a random set of points
in the plane.
"""
import numpy as np
import pylab as pl

from scikits.learn import glm

# this is our test set, it's just a straight line with some
# gaussian noise
xmin, xmax = -5, 5
nsamples = 100
X = [[i] for i in np.linspace(xmin, xmax, nsamples)]
Y = 2 + 0.5 * np.linspace(xmin, xmax, nsamples) +  np.random.randn(nsamples, 1).ravel()

# run the classifier
clf = glm.LinearRegression()
clf.fit(X, Y)

# and plot the result
pl.scatter(X, Y, color='black')
pl.plot(X, clf.predict(X), color='blue', linewidth=3)
pl.show()

