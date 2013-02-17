"""
======================================
Decision Tree Regression with AdaBoost
======================================

A decision tree is boosted using the AdaBoost.R2 [1] algorithm on a 1D
sinusoidal dataset with a small amount of Gaussian noise.
299 boosts (300 decision trees) is compared with a single decision tree
regressor. As the number of boosts is increased the regressor can fit more
detail.

.. [1] H. Drucker, "Improving Regressors using Boosting Techniques", 1997.

"""
print(__doc__)

import numpy as np

# Create a the dataset
rng = np.random.RandomState(1)
X = np.linspace(0, 6, 100)[:, np.newaxis]
y = np.sin(X).ravel() + np.sin(6 * X).ravel() + rng.normal(0, 0.1, X.shape[0])

# Fit regression model
from sklearn.tree import DecisionTreeRegressor
from sklearn.ensemble import AdaBoostRegressor

clf_1 = DecisionTreeRegressor(max_depth=4)

clf_2 = AdaBoostRegressor(DecisionTreeRegressor(max_depth=4),
                          n_estimators=300, random_state=rng)

clf_1.fit(X, y)
clf_2.fit(X, y)

# Predict
y_1 = clf_1.predict(X)
y_2 = clf_2.predict(X)

# Plot the results
import pylab as pl

pl.figure()
pl.scatter(X, y, c="k", label="training samples")
pl.plot(X, y_1, c="g", label="n_estimators=1", linewidth=2)
pl.plot(X, y_2, c="r", label="n_estimators=300", linewidth=2)
pl.xlabel("data")
pl.ylabel("target")
pl.title("Boosted Decision Tree Regression")
pl.legend()
pl.show()
