#!/usr/bin/python
# -*- coding: utf-8 -*-

"""
========================================================================
Gaussian Processes regression: goodness-of-fit on the 'diabetes' dataset
========================================================================

This example consists in fitting a Gaussian Process model onto the diabetes
dataset.
WARNING: This is quite time consuming (about 2 minutes overall runtime).

The corelation parameters are determined by means of maximum likelihood
estimation. An anisotropic squared exponential correlation model with a
constant regression model are assumed. We also used a nugget = 1e-2 in order to
account for the (strong) noise in the targets.

The figure is a goodness-of-fit plot obtained using leave-one-out predictions
of the Gaussian Process model. Based on these predictions, we compute a
leave-one-out of the coefficient of determination (Q2).
"""

# Author: Vincent Dubourg <vincent.dubourg@gmail.com>
# License: BSD style

from scikits.learn import datasets
from scikits.learn.gaussian_process import GaussianProcess
from matplotlib import pyplot as pl

# Print the docstring
print __doc__

# Load the dataset from scikits' data sets
diabetes = datasets.load_diabetes()
X, y = diabetes['data'], diabetes['target']

# Instanciate a GP model
gp = GaussianProcess(regr='constant', corr='absolute_exponential',
                     theta0=[1e-4] * 10, thetaL=[1e-12] * 10,
                     thetaU=[1e-2] * 10, nugget=1e-2, optimizer='Welch',
                     verbose=True)

# Fit the GP model to the data
gp.fit(X, y)

# Estimate the leave-one-out coefficient of determination score
Q2, y_pred = gp.score(return_predictions=True)

# Goodness-of-fit plot
pl.figure()
pl.title('Goodness-of-fit plot (Q2 = %1.2e)' % Q2)
pl.plot(y, y_pred, 'r.', label='Leave-one-out')
pl.plot(y, gp.predict(X), 'k.', label='Whole dataset (nugget=1e-2)')
pl.plot([y.min(), y.max()], [y.min(), y.max()], 'k--')
pl.xlabel('Observations')
pl.ylabel('Predictions')
pl.legend(loc='upper left')
pl.axis('tight')
pl.show()
