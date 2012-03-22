"""
======================
Ordinary Least Squares
======================

Simple Ordinary Least Squares example, we draw the linear least
squares solution for a random set of points in the plane.
"""
print __doc__

import pylab as pl

from sklearn import linear_model
from sklearn.datasets.samples_generator import make_regression

# this is our test set, it's just a straight line with some
# gaussian noise
X, Y = make_regression(n_samples=100, n_features=1, n_informative=1,\
                        random_state=0, noise=35)

# run the classifier
clf = linear_model.LinearRegression()
clf.fit(X, Y)

# and plot the result
pl.scatter(X, Y, color='black')
pl.plot(X, clf.predict(X), color='blue', linewidth=3)
pl.xticks(())
pl.yticks(())
pl.show()
