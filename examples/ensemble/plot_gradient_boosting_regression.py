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
from sklearn.metrics import mean_squared_error

###############################################################################
# Load data
boston = datasets.load_boston()
X, y = shuffle(boston.data, boston.target, random_state=13)
X = X.astype(np.float32)
offset = int(X.shape[0] * 0.9)
X_train, y_train = X[:offset], y[:offset]
X_test, y_test = X[offset:], y[offset:]

################################################################################
# Fit regression model
params = {'n_estimators': 500, 'max_depth': 4, 'min_samples_split': 1,
          'learn_rate': 0.01, 'loss': 'ls'}
clf = ensemble.GradientBoostingRegressor(**params)

clf.fit(X_train, y_train)
mse = mean_squared_error(y_test, clf.predict(X_test))
print("MSE: %.4f" % mse)

################################################################################
# Plot training deviance

# compute test set deviance
y_pred = clf.init.predict(X_test)
test_deviance = np.zeros((params['n_estimators'],), dtype=np.float64)
for i, tree in enumerate(clf.estimators_):
    y_pred += clf.learn_rate * tree.predict(X_test).ravel()
    test_deviance[i] = clf.loss_(y_test, y_pred)

pl.figure()  #figsize=(12, 6))
pl.subplot(1, 2, 1)
pl.title('Deviance')
pl.plot(np.arange(params['n_estimators']) + 1, clf.train_deviance, 'b-',
        label='Training Set Deviance')
pl.plot(np.arange(params['n_estimators']) + 1, test_deviance, 'r-',
        label='Test Set Deviance')
pl.legend(loc='upper right')
pl.xlabel('Boosting Iterations')
pl.ylabel('Deviance')

################################################################################
# Plot feature importance
feature_importance = clf.feature_importances_
sorted_idx = np.argsort(feature_importance)
pos = np.arange(sorted_idx.shape[0]) + .5
pl.subplot(1, 2, 2)
pl.barh(pos, feature_importance[sorted_idx], align='center')
pl.yticks(pos, boston.feature_names[sorted_idx])
pl.xlabel('Relative Importance')
pl.title('Variable Importance')
pl.show()
