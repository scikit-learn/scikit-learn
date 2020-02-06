"""
============================
Gradient Boosting regression
============================

This example demonstrate Gradient Boosting. Gradient boosting can be used for
regression and classification problems. Here, we will use Diabetes (regression)
and breast cancer (classification) datasets.
We will obtain the results from
:class:`~sklearn.ensemble.GradientBoostingRegressor` with least squares loss
and 500 regression trees of depth 4.

"""
print(__doc__)

# Author: Peter Prettenhofer <peter.prettenhofer@gmail.com>
#         Maria Telenczuk <https://github.com/maikia>
#
# License: BSD 3 clause

import numpy as np
import matplotlib.pyplot as plt

from sklearn import ensemble
from sklearn import datasets
from sklearn.utils import shuffle
from sklearn.metrics import mean_squared_error

##############################################################################
# Load the data
# -------------------------------------
#
# First we need to load the data. We will work on the two datasets
# simultaneously.

diabetes = datasets.load_diabetes()
cancer = datasets.load_breast_cancer()


# 
#  
# 
# Center target to avoid gradient boosting init bias: gradient
# boosting
# with the 'recursion' method does not account for the initial estimator
# (here the average target, by default)


Xd, yd = shuffle(diabetes.data, diabetes.target, random_state=13)
X, y = shuffle(cancer.data, cancer.target, random_state=13)
Xd = Xd.astype(np.float32)
X = X.astype(np.float32)
offset = int(X.shape[0] * 0.9)
X_train, y_train = X[:offset], y[:offset]
X_test, y_test = X[offset:], y[offset:]

# #############################################################################
# Fit regression model
params = {'n_estimators': 500, 'max_depth': 4, 'min_samples_split': 2,
          'learning_rate': 0.01, 'loss': 'ls'}
clf = ensemble.GradientBoostingRegressor(**params)

clf.fit(X_train, y_train)
mse = mean_squared_error(y_test, clf.predict(X_test))
print("MSE: %.4f" % mse)

# #############################################################################
# Plot training deviance

# compute test set deviance
test_score = np.zeros((params['n_estimators'],), dtype=np.float64)

for i, y_pred in enumerate(clf.staged_predict(X_test)):
    test_score[i] = clf.loss_(y_test, y_pred)

plt.figure(figsize=(12, 6))
plt.subplot(1, 2, 1)
plt.title('Deviance')
plt.plot(np.arange(params['n_estimators']) + 1, clf.train_score_, 'b-',
         label='Training Set Deviance')
plt.plot(np.arange(params['n_estimators']) + 1, test_score, 'r-',
         label='Test Set Deviance')
plt.legend(loc='upper right')
plt.xlabel('Boosting Iterations')
plt.ylabel('Deviance')

# #############################################################################
# Plot impurity-based feature importance
#
# Warning: impurity-based feature importances can be misleading for
# high cardinality features (many unique values). See
# :func:`sklearn.inspection.permutation_importance` as an alternative.

feature_importance = clf.feature_importances_
# make importances relative to max importance
feature_importance = 100.0 * (feature_importance / feature_importance.max())
sorted_idx = np.argsort(feature_importance)
pos = np.arange(sorted_idx.shape[0]) + .5
plt.subplot(1, 2, 2)
plt.barh(pos, feature_importance[sorted_idx], align='center')
#plt.yticks(pos, np.array(diabetes.feature_names)[sorted_idx])
plt.yticks(pos, np.array(cancer.feature_names)[sorted_idx])
plt.xlabel('Relative Importance')
plt.title('Variable Importance')
plt.show()
