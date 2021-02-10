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
from tqdm import tqdm

np.random.seed(1)


def f(x):
    """The function to predict."""
    return x * np.sin(x)

#----------------------------------------------------------------------
#  First the noiseless case
X = np.atleast_2d(np.random.uniform(0, 10.0, size=100)).T
X = X.astype(np.float32)

# Observations
y = f(X).ravel()

dy = 1.5 + 1.0 * np.random.random(y.shape)
noise = np.random.normal(0, dy)
y += noise
y = y.astype(np.float32)

# Split into train, test datasets.
X_train, X_test, y_train, y_test = train_test_split(X, y)

# Mesh the input space for evaluations of the real function, the prediction and
# its MSE.
xx = np.atleast_2d(np.linspace(0, 10, 1000)).T
xx = xx.astype(np.float32)

# %%
# Let's fit a model based trained with MSE loss.
gbr_mse = GradientBoostingRegressor(loss='ls', n_estimators=250, max_depth=3,
                                    learning_rate=.1, min_samples_leaf=9,
                                    min_samples_split=9)

# %%
# Let's fit this models and three others trained with
# the quantile loss and alpha=0.05, 0.5, 0.95.
# The models obtained for alpha=0.05 and alpha=0.95
# produce a 95% confidence interval. The model trained
# with alpha=0.5 produce the median regression:
# there are the same number of targets above and below
# the predicted values.
gbrs = {
    "q%1.2f" % alpha:
        GradientBoostingRegressor(loss='quantile', alpha=alpha,
                                  n_estimators=250, max_depth=3,
                                  learning_rate=.1, min_samples_leaf=9,
                                  min_samples_split=9)
    for alpha in [0.05, 0.5, 0.95]
}
gbrs['MSE'] = gbr_mse

for alpha, gbr in tqdm(gbrs.items()):
    gbr.fit(X_train, y_train)

# %%
# Let's measure the models given :func:`mean_square_error` and
# :func:`pinball_loss` metrics on the training datasets.
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
# The number of the diagonal is the metric evaluated on a
# model trained to minimize this same metric. It is expected
# to be the minimum value on this column.

# %%
# We do the same on the test set.
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
# Let's finally plot the function, the prediction and the 
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
