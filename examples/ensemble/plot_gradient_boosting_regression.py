"""
============================
Gradient Boosting regression
============================

Demonstrate Gradient Boosting on the boston housing dataset.

This example fits a Gradient Boosting model with least squares loss and
100 regression trees of depth 4. 
"""
print __doc__

# Author: Peter Prettenhofer <peter.prettenhofer@gmail.com>
#
# License: BSD

import numpy as np
import pylab as pl
from sklearn import ensemble
from sklearn import datasets
from sklearn.utils import shuffle
from sklearn import cross_validation

###############################################################################
# Load data
boston = datasets.load_boston()
X, y = shuffle(boston.data, boston.target, random_state=13)
offset = int(X.shape[0] * 0.9)
X_train, y_train = X[:offset], y[:offset]
X_test, y_test = X[offset:], y[offset:]

################################################################################
# Fit regression model
params = {'n_iter': 100, 'max_depth': 4, 'min_split': 1, 'learn_rate': 0.1,
          'loss': 'ls'}
clf = ensemble.GradientBoostingRegressor(**params)

clf.fit(X_train, y_train)
mse = np.mean((clf.predict(X_test) - y_test) ** 2.0)

pl.figure(figsize=(12, 6))
pl.subplot(1,2,1)
pl.title('Deviance')
pl.plot(np.arange(params['n_iter']) + 1, clf.train_deviance, "r-")
pl.xlabel('Boosting Iterations')
pl.ylabel('Training Set Deviance')

################################################################################
# get variable importance
variable_importance = clf.variable_importance
sorted_idx = np.argsort(variable_importance)
pos = np.arange(sorted_idx.shape[0]) + .5
pl.subplot(1,2,2)
pl.barh(pos, variable_importance[sorted_idx], align='center')
pl.yticks(pos, boston.feature_names[sorted_idx])
pl.xlabel('Relative importance')
pl.title('Relative importance of the features for the Boston housing data.')
pl.show()
