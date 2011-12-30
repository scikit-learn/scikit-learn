"""
===================================================================
Boosted Decision Tree Regression
===================================================================

1D regression with boosted :ref:`decision trees <tree>`: the decision tree is
used to fit a sine curve with addition noisy observation. As a result, it
learns local linear regressions approximating the sine curve.

We can see that if the maximum number of boosts (controlled by the
`n_estimators` parameter) is set too high, the ensemble learns too fine
details of the training data and learn from the noise, i.e. they overfit.
"""
print __doc__

import numpy as np

# Create a random dataset
rng = np.random.RandomState(1)
X = np.sort(5 * rng.rand(80, 1), axis=0)
y = np.sin(X).ravel()
y[::5] += 3 * (0.5 - rng.rand(16))

# Fit regression model
from sklearn.tree import DecisionTreeRegressor
from sklearn.ensemble import AdaBoostRegressor

clf_1 = AdaBoostRegressor(DecisionTreeRegressor(max_depth=3), n_estimators=1)
clf_2 = AdaBoostRegressor(DecisionTreeRegressor(max_depth=3), n_estimators=10)

clf_1.fit(X, y)
clf_2.fit(X, y)

# Predict
X_test = np.arange(0.0, 5.0, 0.01)[:, np.newaxis]
y_1 = clf_1.predict(X_test)
y_2 = clf_2.predict(X_test)

# Plot the results
import pylab as pl

pl.figure()
pl.scatter(X, y, c="k", label="data")
pl.plot(X_test, y_1, c="g", label="n_estimators=1", linewidth=2)
pl.plot(X_test, y_2, c="r", label="n_estimators=10", linewidth=2)
pl.xlabel("data")
pl.ylabel("target")
pl.title("Boosted Decision Tree Regression")
pl.legend()
pl.show()
