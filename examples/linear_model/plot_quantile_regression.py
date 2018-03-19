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

from sklearn.linear_model import QuantileRegressor

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
