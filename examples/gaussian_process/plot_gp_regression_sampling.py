#!/usr/bin/python
# -*- coding: utf-8 -*-

"""
=========================================================
Gaussian Processes regression: sampling from the GP
=========================================================

A simple one-dimensional regression exercise with a cubic correlation
model whose parameters are estimated using the maximum likelihood principle.

The figure illustrates the stochastic nature of the Gaussian Processes,
it shows the estimated mean function, which is also the predictor,
a set of functions sampled from the GP stochastic process,
the predicted covariance matrix the predicted 95% variance range.
"""
print __doc__

# Author: Demian Wassermann <demian@bwh.harvard.edu>
# License: BSD style

import numpy as np
from scikits.learn.gaussian_process import GaussianProcess
from matplotlib import pyplot as pl


def f(x):
    """The function to predict."""
    return x * np.sin(x)

# The design of experiments
X = np.atleast_2d([1., 3., 5., 6., 7., 8.]).T

# Observations
y = f(X).ravel()

# Mesh the input space for evaluations of the real function, the prediction and
# its MSE
x = np.atleast_2d(np.linspace(0, 10, 1000)).T

# Instanciate a Gaussian Process model
gp = GaussianProcess(corr='cubic', theta0=1e-2, thetaL=1e-4, thetaU=1e-1, \
                     random_start=100)

# Fit to data using Maximum Likelihood Estimation of the parameters
gp.fit(X, y)

# Make the prediction on the meshed x-axis (ask for MSE as well)
y_pred = gp.predict(x)
covariance_pred = gp.predict_covariance_matrix(x)
sigma_pred = np.sqrt(covariance_pred.diagonal())
y_samples = gp.sample(x, size=10)

# Plot the function, the prediction and the 95% predictive variance
fig = pl.figure()

pl.hold(True)

pl.plot(X, y, 'r.', markersize=15, label=u'Observations')
pl.plot(x, y_pred, 'b-', lw=2, label=u'Prediction')
pl.plot(x, f(x), 'r:', lw=3, label=u'$f(x) = x\,\sin(x)$')

pl.fill(np.concatenate([x, x[::-1]]), \
        np.concatenate([y_pred - 1.9600 * sigma_pred,
                       (y_pred + 1.9600 * sigma_pred)[::-1]]), \
        alpha=.5, fc='b', ec='None', label='95\% predicted variance')

#Plot the samples
pl.plot(x, y_samples[0], '--', c=pl.cm.spring(0.), label=u'Sampled function')
for i, y_sample in enumerate(y_samples[1:]):
    pl.plot(x, y_sample, '--', c=pl.cm.spring((i + 1.) / len(y_samples)))


pl.xlabel('$x$')
pl.ylabel('$f(x)$')
pl.ylim(-10, 20)
pl.legend(loc='upper left')
pl.show()
