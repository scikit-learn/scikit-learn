#!/usr/bin/python
# -*- coding: utf-8 -*-

r"""
==============================================================
Comparison of predictive distributions of different regressors
==============================================================

A simple one-dimensional regression problem adressed by three different
regressors:

1. A Gaussian Process
2. A Random Forest
3. A Bagging-based Regressor

The regressors are fitted based on noisy observations where the magnitude of
the noise at the different training point varies and is known. Plotted are
both the mean and the pointwise 95% confidence interval of the predictions.

This example is based on the example gaussian_process/plot_gp_regression.py.

"""
print(__doc__)

# Author: Jan Hendrik Metzen <jhm@informatik.uni-bremen.de>
# Licence: BSD 3 clause

import numpy as np
from sklearn.gaussian_process import GaussianProcess
from sklearn.ensemble import RandomForestRegressor, BaggingRegressor
from matplotlib import pyplot as pl

np.random.seed(1)


def f(x):
    """The function to predict."""
    return x * np.sin(x)

X = np.linspace(0.1, 9.9, 20)
X = np.atleast_2d(X).T

# Observations and noise
y = f(X).ravel()
dy = 0.5 + 1.0 * np.random.random(y.shape)
noise = np.random.normal(0, dy)
y += noise

# Mesh the input space for evaluations of the real function, the prediction and
# its standard deviation
x = np.atleast_2d(np.linspace(0, 10, 1000)).T

regrs = {"gaussian_process": GaussianProcess(corr='squared_exponential',
                                             theta0=1e-1, thetaL=1e-3,
                                             thetaU=1, nugget=(dy / y) ** 2,
                                             random_start=100),
         "random_forest": RandomForestRegressor(n_estimators=250),
         "bagging": BaggingRegressor(n_estimators=250)}


# Plot predictive distributions of different regressors
fig = pl.figure()
colors = {"gaussian_process": 'b', "random_forest": 'g', "bagging": 'c'}
for name, regr in regrs.items():
    regr.fit(X, y)

    # Make the prediction on the meshed x-axis (ask for standard deviation
    # as well)
    y_pred, sigma = regr.predict(x, with_std=True)

    # Plot 95% confidence interval based on the predictive standard deviation
    pl.plot(x, y_pred, colors[name], label=name)
    pl.fill(np.concatenate([x, x[::-1]]),
            np.concatenate([y_pred - 1.9600 * sigma,
                           (y_pred + 1.9600 * sigma)[::-1]]),
            alpha=.3, fc=colors[name], ec='None')

# Plot the function and the observations
pl.plot(x, f(x), 'r:', label=u'$f(x) = x\,\sin(x)$')
pl.errorbar(X.ravel(), y, dy, fmt='r.', markersize=10, label=u'Observations')
pl.xlabel('$x$')
pl.ylabel('$f(x)$')
pl.ylim(-10, 20)
pl.legend(loc='upper left')

pl.show()
