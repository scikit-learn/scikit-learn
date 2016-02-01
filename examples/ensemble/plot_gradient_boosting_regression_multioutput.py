"""
=========================================
Multi-output Gradient Boosting Regression
=========================================

An example to illustrate multi-output regression with gradient boosting.

This example illustrates the use of the
:ref:`multioutput.MultiOutputRegressor <_multiclass>` meta-estimator
to perform multi-output regression with an estimator that does not
natively supoort it. For comparison a regression tree based model is
shown which supports multi-output regression natively.

Using a single underlying feature the model learns both the
x and y coordinate as output.

"""
print(__doc__)

# Author: Tim Head <betatim@gmail.com>
#
# License: BSD 3 clause

import numpy as np
import matplotlib.pyplot as plt
from sklearn.ensemble import GradientBoostingRegressor
from sklearn.multioutput import MultiOutputRegressor
from sklearn.tree import DecisionTreeRegressor

# Create a random dataset
rng = np.random.RandomState(1)
X = np.sort(200 * rng.rand(200, 1) - 100, axis=0)
y = np.array([np.pi * np.sin(X).ravel(), np.pi * np.cos(X).ravel()]).T
y += (0.5 - rng.rand(*y.shape))

# Fit gradient boosted trees and a single tree
regr_gbr = MultiOutputRegressor(GradientBoostingRegressor(max_depth=6, random_state=0))
regr_gbr.fit(X, y)

regr_tree = DecisionTreeRegressor(max_depth=6, random_state=2)
regr_tree.fit(X, y)

# Predict
X_test = np.arange(-100.0, 100.0, 0.01)[:, np.newaxis]
y_gbr = regr_gbr.predict(X_test)
y_tree = regr_tree.predict(X_test)

# Plot the results
plt.figure()
s = 50
plt.scatter(y[:, 0], y[:, 1], c="navy", s=s, marker="^", label="Data")
plt.scatter(y_gbr[:, 0], y_gbr[:, 1], c="cornflowerblue", s=s, label="GBR")
plt.scatter(y_tree[:, 0], y_tree[:, 1], c="c", s=s, marker="s", label="Tree")
plt.xlim([-6, 6])
plt.ylim([-6, 6])
plt.xlabel("target 1")
plt.ylabel("target 2")
plt.title("Multi-output Gradient Boosting Regression")
plt.legend()
plt.show()
