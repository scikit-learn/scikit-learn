#!/usr/bin/python
# -*- coding: utf-8 -*-

r"""
==============================================================
Comparison of predictive distributions of different regressors
==============================================================

A simple one-dimensional, noisy regression problem adressed by three different
regressors:

1. A Gaussian Process
2. A Random Forest
3. A Bagging-based Regressor

The regressors are fitted based on noisy observations where the magnitude of
the noise at the different training point is constant and known. Plotted are
both the mean and the pointwise 95% confidence interval of the predictions.
The mean predictions are evaluated on noise-less test data using the mean-
squared-error. The mean log probabilities of the noise-less test data are used
to evaluate the predictive distributions (a normal distribution with the
predicted mean and standard deviation) of the three regressors.

This example is based on the example gaussian_process/plot_gp_regression.py.
"""
print(__doc__)

# Author: Jan Hendrik Metzen <jhm@informatik.uni-bremen.de>
# Licence: BSD 3 clause

import numpy as np
from scipy.stats import norm
from sklearn.gaussian_process import GaussianProcess
from sklearn.ensemble import RandomForestRegressor, BaggingRegressor
from sklearn.metrics import mean_squared_error
from matplotlib import pyplot as pl

np.random.seed(1)


def f(x):
    """The function to predict."""
    return x * np.sin(x)

X = np.linspace(0.1, 9.9, 20)
X = np.atleast_2d(X).T

# Observations and noise
y = f(X).ravel()
dy = np.ones_like(y)
noise = np.random.normal(0, dy)
y += noise

# Mesh the input space for evaluations of the real function, the prediction and
# its standard deviation
x = np.atleast_2d(np.linspace(0, 10, 1000)).T

regrs = {"Gaussian Process": GaussianProcess(corr='squared_exponential',
                                             theta0=1e-1, thetaL=1e-3,
                                             thetaU=1, nugget=(dy / y) ** 2,
                                             random_start=100),
         "Random Forest": RandomForestRegressor(n_estimators=250),
         "Bagging": BaggingRegressor(n_estimators=250)}


# Plot predictive distributions of different regressors
fig = pl.figure()
# Plot the function and the observations
pl.plot(x, f(x), 'r', label=u'$f(x) = x\,\sin(x)$')
pl.fill(np.concatenate([x, x[::-1]]),
        np.concatenate([f(x) - 1.9600, (f(x) + 1.9600)[::-1]]),
        alpha=.3, fc='r', ec='None')
pl.plot(X.ravel(), y, 'ko', zorder=5, label=u'Observations')
# Plot predictive distibutions of GP and Bagging
colors = {"Gaussian Process": 'b', "Bagging": 'g'}
mse = {}
log_pdf_loss = {}
for name, regr in regrs.items():
    regr.fit(X, y)

    # Make the prediction on the meshed x-axis (ask for standard deviation
    # as well)
    y_pred, sigma = regr.predict(x, with_std=True)

    # Compute mean-squared error and log predictive loss
    mse[name] = mean_squared_error(f(x), y_pred)
    log_pdf_loss[name] = \
        norm(y_pred, sigma).logpdf(f(x)).mean()

    if name == "Random Forest":  # Skip because RF is very similar to Bagging
        continue

    # Plot 95% confidence interval based on the predictive standard deviation
    pl.plot(x, y_pred, colors[name], label=name)
    pl.fill(np.concatenate([x, x[::-1]]),
            np.concatenate([y_pred - 1.9600 * sigma,
                           (y_pred + 1.9600 * sigma)[::-1]]),
            alpha=.3, fc=colors[name], ec='None')


pl.xlabel('$x$')
pl.ylabel('$f(x)$')
pl.ylim(-10, 20)
pl.legend(loc='upper left')

print "Mean-squared error of predictors on 1000 equidistant noise-less test " \
    "datapoints:\n\tRandom Forest: %.2f\n\tBagging: %.2f" \
    "\n\tGaussian Process: %.2f" \
    % (mse["Random Forest"], mse["Bagging"], mse["Gaussian Process"])

print "Mean log-probability of 1000 equidistant noise-less test datapoints " \
    "under the (normal) predictive distribution of the predictors, i.e., " \
    "log N(y_true| y_pred_mean, y_pred_std) [less is better]:"\
    "\n\tRandom Forest: %.2f\n\tBagging: %.2f\n\tGaussian Process: %.2f" \
    % (log_pdf_loss["Random Forest"], log_pdf_loss["Bagging"],
       log_pdf_loss["Gaussian Process"])

print "In summary, the mean predictions of the Gaussian Process are slightly "\
    "better than those of Random Forest and Bagging. The predictive " \
    "distributions (taking into account also the predictive variance) " \
    "of the Gaussian Process are considerably better."

pl.show()
