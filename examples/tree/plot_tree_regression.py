"""
===================================================================
Decision Tree Regression
===================================================================

Toy example of 1D regression.

"""
print __doc__

###############################################################################
# Generate sample data
import numpy as np

X = np.sort(5*np.random.rand(40, 1), axis=0)
y = np.sin(X).ravel()

###############################################################################
# Add noise to targets
y[::5] += 3*(0.5 - np.random.rand(8))

###############################################################################
# Fit regression model
from sklearn.tree import DecisionTreeRegressor

clf_1 = DecisionTreeRegressor(max_depth=2)
clf_2 = DecisionTreeRegressor(max_depth=5)
y_1 = clf_1.fit(X, y).predict(X)
y_2 = clf_2.fit(X, y).predict(X)

###############################################################################
# look at the results
import pylab as pl
pl.scatter(X, y, c='k', label='data')
pl.hold('on')
pl.plot(X, y_1, c='g', label='max_depth=2')
pl.plot(X, y_2, c='r', label='max_depth=5')
pl.xlabel('data')
pl.ylabel('target')
pl.title('Decision Tree Regression')
pl.legend()
pl.show()