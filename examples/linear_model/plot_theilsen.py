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

from __future__ import division, print_function, absolute_import
import numpy as np
import matplotlib.pylab as pl
from sklearn.linear_model import LinearRegression, TheilSen, RANSACRegressor

print(__doc__)

###################
# Theil-Sen / OLS #
###################

# Generate 1d toy problem
np.random.seed(0)
n_samples = 100
# Linear model y = 3*x + N(2, 0.1**2)
x = np.random.randn(n_samples)
w = np.array([3.])
c = np.array([2.])
noise = 0.1 * np.random.randn(n_samples)
y = w * x + c + noise
# Add some outliers
x[42], y[42] = (-2, 4)
x[43], y[43] = (-2.5, 8)
x[53], y[53] = (2.5, 1)
x[60], y[60] = (2.1, 2)
x[72], y[72] = (1.8, -7)
X = x[:, np.newaxis]

# Ordinary Least Squares fit
lstq = LinearRegression().fit(X, y)
line_x = np.array([-3, 3])
y_pred_lstq = lstq.predict(line_x.reshape(2, 1))
# Theil-Sen fit
theilsen = TheilSen().fit(X, y)
y_pred_theilsen = theilsen.predict(line_x.reshape(2, 1))
# Plot
fig, ax = pl.subplots()
pl.scatter(X, y)
pl.plot(line_x, y_pred_lstq, label='OLS')
pl.plot(line_x, y_pred_theilsen, label='Theil-Sen')
ax.set_xlim(line_x)
pl.legend(loc=2)

######################
# Theil-Sen / RANSAC #
######################

# Construct special case
x = np.array([0.16, 0.16, 0.06, 0.87, 0.29,
              0.28, 0.22, 0.11, 0.86, 0.34])
y = np.array([23.46, 29.91, 6.66, 99, 52.55,
              44.9, 34.44, 15.31, 98, 61.34])
X = x[:, np.newaxis]
pred_X = np.array([0., 1.]).reshape(2, 1)
pl.figure()
pl.ylim(0, 100)
pl.plot(x, y, 'bo')
# RANSAC fit
rs = RANSACRegressor(random_state=42).fit(X, y)
pred_y = rs.predict(pred_X)
pl.plot(pred_X, pred_y, label="RANSAC")
# Theil-Sen fit
ts = TheilSen(random_state=42).fit(X, y)
pred_y = ts.predict(pred_X)
pl.plot(pred_X, pred_y, label="Theil-Sen")
pl.legend(loc=2)
pl.show()
