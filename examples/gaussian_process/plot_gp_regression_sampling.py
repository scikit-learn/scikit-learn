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
the predicted covariance matrix and the estimated error as the
squared diagonal of the predicted covariance matrix in the form of a
pointwise 95% confidence interval.
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
gp = GaussianProcess(corr='squared_exponential', theta0=1e-1, thetaL=1e-4, thetaU=5, \
                     random_start=100)

# Fit to data using Maximum Likelihood Estimation of the parameters
gp.fit(X, y)

# Make the prediction on the meshed x-axis (ask for MSE as well)
y_pred = gp.predict(x)
covariance_pred = gp.predict_covariance_matrix(x)
sigma_pred = np.sqrt(covariance_pred.diagonal())
y_samples = gp.sample(x, size=10)

# Plot the function, the prediction and the 95% confidence interval based on
# the MSE
fig = pl.figure()
pl.subplot(2,1,1)

pl.hold(True)

#Plot the samples
for y_sample in y_samples:
    pl.plot(x, y_sample,'--')

pl.plot(X, y, 'r.', markersize=10, label=u'Observations')
pl.plot(x, y_pred, 'b-', lw=3, label=u'Prediction')
pl.plot(x, f(x), 'r:', lw=3, label=u'$f(x) = x\,\sin(x)$')

#pl.fill(np.concatenate([x, x[::-1]]), \
#        np.concatenate([y_pred - 1.9600 * sigma,
#                       (y_pred + 1.9600 * sigma)[::-1]]), \
#        alpha=.5, fc='b', ec='None', label='95% confidence interval')
pl.fill(np.concatenate([x, x[::-1]]), \
        np.concatenate([y_pred - 1.9600 * sigma_pred,
                       (y_pred + 1.9600 * sigma_pred)[::-1]]), \
        alpha=.5, fc='g', ec='None', label='95% predicted variance')

pl.xlabel('$x$')
pl.ylabel('$f(x)$')
pl.ylim(-10, 20)
pl.legend(loc='upper left')

pl.subplot(2,1,2)
pl.imshow(covariance_pred)

pl.show()
