"""
============================
Gradient Boosting regression
============================

This example demonstrate Gradient Boosting which will produce a prediction
model from ensemble of weak prediction models. Gradient boosting can be used
for regression and classification problems. Here, we will use Diabetes
(regression) and breast cancer (classification) datasets.

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
from sklearn.model_selection import train_test_split

##############################################################################
# Load the data
# -------------------------------------
#
# First we need to load the data. We will work on both datasets
# simultaneously. We set random state to be consistent with the result.

diabetes = datasets.load_diabetes()
cancer = datasets.load_breast_cancer()

Xd, yd = diabetes.data, diabetes.target
Xc, yc = cancer.data, cancer.target

##############################################################################
# Data preprocessing
# -------------------------------------
#
# Next, we will split our datasets to 90% for training and leave the rest for
# testing. We will also prepare the parameters we want to use to fit our
# regression model. You can play with those parameters to see how the
# results change:
#
# Here:
# n_estimators : is the number of boosting stages which will be performed.
#     Later, we will plot and see how the deviance changes with those boosting
#     operations.
# max_depth : this limits the number of nodes in the tree. The best value
#     depends on the interaction of the input variables.
# min_samples_split : is the minimum number of samples required to split an
#     internal node.
# learning_rate: tells how much the contribution of each tree will shrink
# loss: here, we decided to use least squeares as a loss function, however
#     there are many other options (check
#     :class:`~sklearn.ensemble.GradientBoostingRegressor` to see what are
#     other possibilities)
# In this example we will use the same parameters for both datasets even if
# they might not be the most optimal

Xd_train, Xd_test, yd_train, yd_test = train_test_split(Xd, yd,
                                                        test_size=0.1,
                                                        random_state=13)

Xc_train, Xc_test, yc_train, yc_test = train_test_split(Xc, yc,
                                                        test_size=0.1,
                                                        random_state=13)

params = {'n_estimators': 500,
          'max_depth': 4,
          'min_samples_split': 5,
          'learning_rate': 0.01,
          'loss': 'ls'}

##############################################################################
# Fit regression model
# -------------------------------------
#
# Now we will initiate the gradient boosting regressors and fit it with our
# training data. Let's also look and the mean squared error on the test data.

clf_d = ensemble.GradientBoostingRegressor(**params)
clf_c = ensemble.GradientBoostingRegressor(**params)

clf_d.fit(Xd_train, yd_train)
clf_c.fit(Xc_train, yc_train)

mse_d = mean_squared_error(yd_test, clf_d.predict(Xd_test))
mse_c = mean_squared_error(yc_test, clf_c.predict(Xc_test))
print("The mean squared error (MSE) on the")
print("diabetes dataset: {:.4f}".format(mse_d))
print("cancer dataset: {:.4f}".format(mse_c))

##############################################################################
# Plot training deviance
# -------------------------------------
#
# Finally, we will visualize the results. To do that we will first compute the
# test set deviance and then plot it. We are going to plot the results for both
# datasets on a single figure.


test_score_c = np.zeros((params['n_estimators'],), dtype=np.float64)
test_score_d = np.zeros((params['n_estimators'],), dtype=np.float64)

for i, yc_pred in enumerate(clf_c.staged_predict(Xc_test)):
    test_score_c[i] = clf_c.loss_(yc_test, yc_pred)

for i, yd_pred in enumerate(clf_d.staged_predict(Xd_test)):
    test_score_d[i] = clf_d.loss_(yd_test, yd_pred)

plt.figure()

plt.subplot(1, 2, 1)
plt.title('Diabetes dataset')
plt.plot(np.arange(params['n_estimators']) + 1, clf_d.train_score_, 'b-',
         label='Training Set Deviance')
plt.plot(np.arange(params['n_estimators']) + 1, test_score_d, 'r-',
         label='Test Set Deviance')
plt.xlabel('Boosting Iterations')
plt.ylabel('Deviance')

plt.subplot(1, 2, 2)
plt.title('Cancer dataset')
plt.plot(np.arange(params['n_estimators']) + 1, clf_c.train_score_, 'b-',
         label='Training Set Deviance')
plt.plot(np.arange(params['n_estimators']) + 1, test_score_c, 'r-',
         label='Test Set Deviance')
plt.legend(loc='upper right')
plt.xlabel('Boosting Iterations')


# #############################################################################
# Plot impurity-based feature importance
#
# Warning: impurity-based feature importances can be misleading for
# high cardinality features (many unique values). See
# :func:`sklearn.inspection.permutation_importance` as an alternative.

plt.figure()
feature_importance = clf_c.feature_importances_
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
