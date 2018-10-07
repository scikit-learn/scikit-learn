"""
==============
Quantile regression
==============

Plot the prediction of different conditional quantiles.

The left figure shows the case when error distribution is normal,
but variance is not constant.

The right figure shows example of an asymmetric error distribution
(namely, Pareto).
"""
from __future__ import division
print(__doc__)

import numpy as np
import matplotlib.pyplot as plt

from sklearn.linear_model import QuantileRegressor, LinearRegression
from sklearn.metrics import mean_absolute_error, mean_squared_error
from sklearn.model_selection import cross_val_score

fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(15, 5))

rng = np.random.RandomState(42)
x = np.linspace(0, 10, 100)
X = x[:, np.newaxis]
y = 20 + x*2 + rng.normal(loc=0, scale=0.5+0.5*x, size=x.shape[0])
ax1.scatter(x, y)

quantiles = [0.05, 0.5, 0.95]
for quantile in quantiles:
    qr = QuantileRegressor(quantile=quantile, max_iter=10000, alpha=0)
    qr.fit(X, y)
    ax1.plot([0, 10], qr.predict([[0], [10]]))
ax1.xlabel('x')
ax1.ylabel('y')
ax1.title('Quantiles of normal residuals with non-constant variance')
ax1.legend(quantiles)

y = 20 + x * 0.5 + rng.pareto(10, size=x.shape[0])*10
ax2.scatter(x, y)

for quantile in quantiles:
    qr = QuantileRegressor(quantile=quantile, max_iter=10000, alpha=0)
    qr.fit(X, y)
    ax2.plot([0, 10], qr.predict([[0], [10]]))
ax2.xlabel('x')
ax2.ylabel('y')
ax2.title('Quantiles of asymmetrically distributed residuals')
ax2.legend(quantiles)

plt.show()

#########################################################################
#
# The second part of the code shows that LinearRegression minimizes RMSE,
# while QuantileRegressor minimizes MAE, and both do their own job well.

models = [LinearRegression(), QuantileRegressor(alpha=0, max_iter=10000)]
names = ['OLS', 'Quantile']

print('# In-sample performance')
for model_name, model in zip(names, models):
    print(model_name + ':')
    model.fit(X, y)
    mae = mean_absolute_error(model.predict(X), y)
    rmse = np.sqrt(mean_squared_error(model.predict(X), y))
    print('MAE={:.4}  RMSE={:.4}'.format(mae, rmse))
print('\n# Cross-validated performance')
for model_name, model in zip(names, models):
    print(model_name + ':')
    mae = -cross_val_score(model, X, y,
                           scoring='neg_mean_absolute_error').mean()
    rmse = np.sqrt(-cross_val_score(model, X, y,
                                    scoring='neg_mean_squared_error').mean())
    print('MAE={:.4}  RMSE={:.4}'.format(mae, rmse))
