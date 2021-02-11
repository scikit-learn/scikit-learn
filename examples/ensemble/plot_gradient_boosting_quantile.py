"""
=====================================================
Prediction Intervals for Gradient Boosting Regression
=====================================================

This example shows how quantile regression can be used
to create prediction intervals.
"""

import numpy as np
import matplotlib.pyplot as plt

from sklearn.ensemble import GradientBoostingRegressor
from sklearn.metrics import pinball_loss, mean_squared_error
from sklearn.model_selection import train_test_split

from pandas import DataFrame

np.random.seed(1)

# %%
# The function to predict.

def f(x):
    return x * np.sin(x)


# %%
# Draw a set of points `(x, y=f(x))`.
X = np.atleast_2d(np.random.uniform(0, 10.0, size=200)).T
X = X.astype(np.float32)
y = f(X).ravel()

# %%
# A gaussian noise is added to get `(x, y=f(x) + \epsilon)`.
dy = 1.5 + 1.0 * np.random.random(y.shape)
noise = np.random.normal(0, dy)
y += noise
y = y.astype(np.float32)

# %%
# (X, y) is split into train, test datasets.
X_train, X_test, y_train, y_test = train_test_split(X, y)

# %%
# Create a set is a deterministic set of points in [0, 10].
xx = np.atleast_2d(np.linspace(0, 10, 1000)).T
xx = xx.astype(np.float32)

# %%
# Instantiate a model based trained with MSE loss.
gbr_mse = GradientBoostingRegressor(loss='ls', n_estimators=250, max_depth=2,
                                    learning_rate=.05, min_samples_leaf=9,
                                    min_samples_split=9)

# %%
# Fit this model and three others trained with
# the quantile loss and alpha=0.05, 0.5, 0.95.
# The models obtained for alpha=0.05 and alpha=0.95
# produce a 95% confidence interval. The model trained
# with alpha=0.5 produces a regression of the median:
# there are the same number of targets above and below
# the predicted values.
gbrs = {
    "q%1.2f" % alpha:
        GradientBoostingRegressor(loss='quantile', alpha=alpha,
                                  n_estimators=250, max_depth=2,
                                  learning_rate=.05, min_samples_leaf=9,
                                  min_samples_split=9)
    for alpha in [0.05, 0.5, 0.95]
}
gbrs['MSE'] = gbr_mse

for alpha, gbr in gbrs.items():
    gbr.fit(X_train, y_train)

# %%
# Plot the function, the prediction and the
# 95% confidence interval based on the MSE.
y_pred = gbrs['MSE'].predict(xx)
y_lower = gbrs['q0.05'].predict(xx)
y_upper = gbrs['q0.95'].predict(xx)
y_med = gbrs['q0.50'].predict(xx)

fig = plt.figure()
plt.plot(xx, f(xx), 'g:', label=r'$f(x) = x\,\sin(x)$')
plt.plot(X_test, y_test, 'b.', markersize=10, label='Test observations')
plt.plot(xx, y_med, 'r-', label='Prediction Median', color="orange")
plt.plot(xx, y_pred, 'r-', label='Prediction MSE')
plt.plot(xx, y_upper, 'k-')
plt.plot(xx, y_lower, 'k-')
plt.fill(np.concatenate([xx, xx[::-1]]),
         np.concatenate([y_upper, y_lower[::-1]]),
         alpha=.5, fc='b', ec='None', label='95% prediction interval')
plt.xlabel('$x$')
plt.ylabel('$f(x)$')
plt.ylim(-10, 20)
plt.legend(loc='upper left')
plt.show()

# %%
# Measure the models given :func:`mean_square_error` and
# :func:`pinball_loss` metrics on the training dataset.
results = []
for name, gbr in sorted(gbrs.items()):
    metrics = {'model': name}
    y_pred = gbr.predict(X_train)
    for alpha in [0.05, 0.5, 0.95]:
        metrics["pbl=%1.2f" % alpha] = pinball_loss(
            y_train, y_pred, alpha=alpha)
    metrics['MSE'] = mean_squared_error(y_train, y_pred)
    results.append(metrics)
DataFrame(results).set_index('model')

# %%
# One column shows all models evaluated by the same metric.
# The minimum number on a column is obtained when the model
# is trained and measured with the same metric.
# This should be always the case on the training set if the training
# converged. The measures give almost the same results
# when the model is trained with MSE or with the quantile loss and alpha=0.5.
# The random noise added to the data explains that proximity.
# Both trainings would give much more different models if the dataset
# had outliers or if the added noise was not gaussian.
# The quantile loss is less sensitive to outliers.
#
# We then do the same on the test set.
results = []
for name, gbr in sorted(gbrs.items()):
    metrics = {'model': name}
    y_pred = gbr.predict(X_test)
    for alpha in [0.05, 0.5, 0.95]:
        metrics["pbl=%1.2f" % alpha] = pinball_loss(
            y_test, y_pred, alpha=alpha)
    metrics['MSE'] = mean_squared_error(y_test, y_pred)
    results.append(metrics)
DataFrame(results).set_index('model')

# %%
# Errors are higher meaning the model slightly overfitted the data.
# It still shows the minimum of a metric is obtained when
# the model is trained by minimizing this same metric.
