"""
======================================
Decision Tree Regression with AdaBoost
======================================

A decision tree is boosted using the AdaBoost.R2 [1]_ algorithm on a 1D
sinusoidal dataset with a small amount of Gaussian noise.
299 boosts (300 decision trees) is compared with a single decision tree
regressor. As the number of boosts is increased the regressor can fit more
detail.

.. [1] H Drucker, "Improving Regressors using Boosting Techniques", 1997.

"""

# %%
# Creating dummy dataset with a sinusoidal relationship and some gaussian noise
# -----------------------------------------------------------------------------

# Author: Noel Dawe <noel.dawe@gmail.com>
#
# License: BSD 3 clause

import numpy as np

rng = np.random.RandomState(1)
X = np.linspace(0, 6, 100)[:, np.newaxis]
y = np.sin(X).ravel() + np.sin(6 * X).ravel() + rng.normal(0, 0.1, X.shape[0])

# %%
# Training and prediction with AdaBoost and DecisionTree Regressors
# -----------------------------------------------------------------
# The base learner is a DecisionTreeRegressor with `max_depth=4`.
# AdaBoostRegressor will be built with `n_estimators=300` of those base learners.

from sklearn.ensemble import AdaBoostRegressor
from sklearn.tree import DecisionTreeRegressor

regr_1 = DecisionTreeRegressor(max_depth=4)

regr_2 = AdaBoostRegressor(
    DecisionTreeRegressor(max_depth=4), n_estimators=300, random_state=rng
)

regr_1.fit(X, y)
regr_2.fit(X, y)

y_1 = regr_1.predict(X)
y_2 = regr_2.predict(X)

# %%
# Plotting the results
# --------------------

import matplotlib.pyplot as plt
import seaborn as sns

# reshape the colors as 1D array for matplotlib
colors = np.asarray(sns.color_palette("colorblind")).reshape(10, 1, 3)

plt.figure()
plt.scatter(X, y, c=colors[0], label="training samples")
plt.plot(X, y_1, c=colors[1], label="n_estimators=1", linewidth=2)
plt.plot(X, y_2, c=colors[2], label="n_estimators=300", linewidth=2)
plt.xlabel("data")
plt.ylabel("target")
plt.title("Boosted Decision Tree Regression")
plt.legend()
plt.show()
