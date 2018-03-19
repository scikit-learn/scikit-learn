"""
==============
Quantile regression
==============

Plot the prediction of different conditional quantiles.

The left figure shows the case when error distribution is normal,
but variance is not constant.

The right figure shows example of an asymmetric error distribution
(namely, Pareto).

The second part of the code shows that LinearRegression minimizes RMSE,
while QuantileRegressor minimizes MAE, and both do their own job well.
"""
from __future__ import division
print(__doc__)

import numpy as np
import matplotlib.pyplot as plt

from sklearn.linear_model import QuantileRegressor, LinearRegression
from sklearn.metrics import mean_absolute_error, mean_squared_error
from sklearn.model_selection import cross_val_score

plt.figure(figsize=(15, 5))

plt.subplot(121)
np.random.seed(0)
x = np.linspace(0, 10, 100)
X = x[:, np.newaxis]
y = 20 + x*2 + np.random.normal(loc=0, scale=0.5+0.5*x, size=x.shape[0])
plt.scatter(x, y)

quantiles = [0.05, 0.5, 0.95]
for quantile in quantiles:
    qr = QuantileRegressor(quantile=quantile, max_iter=10000, alpha=0)
    qr.fit(X, y)
    plt.plot([0, 10], qr.predict([[0], [10]]))
plt.xlabel('x')
plt.ylabel('y')
plt.title('Quantiles of normal residuals with non-constant variance')
plt.legend(quantiles)

plt.subplot(122)

y = 20 + x * 0.5 + np.random.pareto(10, size=x.shape[0])*10
plt.scatter(x, y)

for quantile in quantiles:
    qr = QuantileRegressor(quantile=quantile, max_iter=10000, alpha=0)
    qr.fit(X, y)
    plt.plot([0, 10], qr.predict([[0], [10]]))
plt.xlabel('x')
plt.ylabel('y')
plt.title('Quantiles of asymmetrically distributed residuals')
plt.legend(quantiles)

plt.show()

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
