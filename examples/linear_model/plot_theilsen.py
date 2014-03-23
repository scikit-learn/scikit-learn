"""
====================
Theil-Sen Regression
====================

Computes a Theil-Sen Regression on a synthetic dataset.

See :ref:`theil_sen_regression` for more information on the regressor.

Compared to the OLS (ordinary least squares) estimator, the Theil-Sen
estimator is robust against outliers. It has a breakdown point of about 29.3%
in case of a simple linear regression which means that it can tolerate
arbitrary corrupted data (outliers) of up to 29.3% in the two-dimensional
case.

The estimation of the model is done by calculating the slopes and intercepts
of a subpopulation of all possible combinations of p subsample points. If an
intercept is fitted, p must be greater than or equal to n_features + 1. The
final slope and intercept is then defined as the spatial median of these
slopes and intercepts.

In certain cases Theil-Sen performs better than :ref:`RANSAC
<ransac_regression>` which is also a robust method. Due to the computational
complexity of Theil-Sen it is recommended to use it only for small problems in
terms of number of samples and features. For larger problems the
``max_subpopulation`` parameter restricts the magnitude of all possible
combinations of p subsample points to a randomly chosen subset and therefore
also limits the runtime. Therefore, Theil-Sen is applicable to larger problems
with the drawback of losing some of its mathematical properties since it works
on a random subset.
"""

# Author: Florian Wilhelm -- <florian.wilhelm@gmail.com>
# License: BSD 3 clause

import time
import numpy as np
import matplotlib.pylab as pl
from sklearn.linear_model import LinearRegression, TheilSen, RANSACRegressor

print(__doc__)

estimators = [('OLS', LinearRegression()),
              ('Theil-Sen', TheilSen()),
              ('RANSAC', RANSACRegressor(random_state=42)), ]

##############################################################################
# Outliers only in the y direction

np.random.seed(0)
n_samples = 200
# Linear model y = 3*x + N(2, 0.1**2)
x = np.random.randn(n_samples)
w = 3.
c = 2.
noise = 0.1 * np.random.randn(n_samples)
y = w * x + c + noise
# 10% outliers
y[-20:] += -20 * x[-20:]
X = x[:, np.newaxis]

pl.plot(x, y, 'k+', mew=2, ms=8)
line_x = np.array([-3, 3])
for name, estimator in estimators:
    t0 = time.time()
    estimator.fit(X, y)
    elapsed_time = time.time() - t0
    y_pred = estimator.predict(line_x.reshape(2, 1))
    pl.plot(line_x, y_pred,
            label='%s (fit time: %.2fs)' % (name, elapsed_time))

pl.axis('tight')
pl.legend(loc=2)


##############################################################################
# Outliers in the X direction

np.random.seed(0)
# Linear model y = 3*x + N(2, 0.1**2)
x = np.random.randn(n_samples)
noise = 0.1 * np.random.randn(n_samples)
y = 3 * x + 2 + noise
# 10% outliers
x[-20:] = 9.9
y[-20:] += 15
X = x[:, np.newaxis]

pl.figure()
pl.plot(x, y, 'k+', mew=2, ms=8)

line_x = np.array([-3, 10])
for name, estimator in estimators:
    t0 = time.time()
    estimator.fit(X, y)
    elapsed_time = time.time() - t0
    y_pred = estimator.predict(line_x.reshape(2, 1))
    pl.plot(line_x, y_pred,
            label='%s (fit time: %.2fs)' % (name, elapsed_time))

pl.axis('tight')
pl.legend(loc=2)
pl.show()
