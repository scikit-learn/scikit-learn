"""
====================
Theil-Sen Regression
====================

Computes a Theil-Sen Regression on a synthetic dataset.

See :ref:`theil_sen_regression` for more information on the regressor.

Compared to the OLS (ordinary least squares) estimator, the Theil-Sen
estimator is robust against outliers. It has a breakdown point of about 29.3%
in case of a simple linear regression which means that it can tolerate
arbitrary corrupted data (outliers) of up to 29.3%.

The estimation of the model is done by calculating the slopes and intercepts
of a subpopulation of all possible combinations of p+1 sample points given
that p is the number of features. The final slope and intercept is then
defined as the spatial median of these slopes and intercepts.
"""

# Author: Florian Wilhelm -- <florian.wilhelm@gmail.com>
# License: BSD 3 clause

from __future__ import division, print_function, absolute_import
import numpy as np
import matplotlib.pylab as plt
from sklearn.linear_model import LinearRegression, TheilSen

print(__doc__)

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
fig, ax = plt.subplots(figsize=(6, 5))
plt.scatter(X, y)
plt.plot(line_x, y_pred_lstq, label='OLS')
plt.plot(line_x, y_pred_theilsen, label='Theil-Sen')
ax.set_xlim(line_x)
plt.legend()
plt.show()