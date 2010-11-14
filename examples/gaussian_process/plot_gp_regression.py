#!/usr/bin/python
# -*- coding: utf-8 -*-

"""
===============================================
Gaussian Processes for Machine Learning example
===============================================

A simple one-dimensional regression exercise with a
cubic correlation model whose parameters are estimated
using the maximum likelihood principle.
"""
print __doc__
# Author: Vincent Dubourg <vincent.dubourg@gmail.com
# License: BSD style

import numpy as np
from scipy import stats
from scikits.learn.gaussian_process import GaussianProcess, corrcubic
from matplotlib import pyplot as pl

# The function to predict
f = lambda x: x*np.sin(x)

# The design of experiments
X = np.array([1., 3., 5., 6., 7., 8.])

# Observations
Y = f(X)

# Mesh the input space for evaluations of the real function, the prediction and its MSE
x = np.linspace(0,10,1000)

# Instanciate a Gaussian Process model
gp = GaussianProcess(corr=corrcubic, theta0=1e-2, thetaL=1e-4, thetaU=1e-1, random_start=100)

# Fit to data using Maximum Likelihood Estimation of the parameters
gp.fit(X, Y)

# Make the prediction on the meshed x-axis (ask for MSE as well)
y, MSE = gp.predict(x, eval_MSE=True)
sigma = np.sqrt(MSE)

# Plot the function, the prediction and the 95% confidence interval based on the MSE
fig = pl.figure()
pl.plot(x, f(x), 'r:', label=u'$f(x) = x\,\sin(x)$')
pl.plot(X, Y, 'r.', markersize=10, label=u'Observations')
pl.plot(x, y, 'b-', label=u'Prediction')
pl.fill(np.concatenate([x, x[::-1]]), np.concatenate([y - 1.9600 * sigma, (y + 1.9600 * sigma)[::-1]]), alpha=.5, fc='b', ec='None', label='95% confidence interval')
pl.xlabel('$x$')
pl.ylabel('$f(x)$')
pl.ylim(-10, 20)
pl.legend(loc='upper left')

pl.show()
