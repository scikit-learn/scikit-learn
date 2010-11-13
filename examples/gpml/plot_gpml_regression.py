#!/usr/bin/python
# -*- coding: utf-8 -*-

"""
===============================================
Gaussian Processes for Machine Learning example
===============================================

A simple one-dimensional regression exercise with
different correlation models and maximum likelihood
estimation of the Gaussian Process Model parameters.
"""
# Author: Vincent Dubourg <vincent.dubourg@gmail.com
# License: BSD style

import numpy as np
from scipy import stats
from scikits.learn.gpml import GaussianProcessModel, correxp1, correxp2, corrlin, corrcubic
from matplotlib import pyplot as pl

# The function to predict
f = lambda x: x*np.sin(x)

# The design of experiments
X = np.array([1., 3., 5., 6., 7., 8.])

# Observations
Y = f(X)

# Mesh the input space for evaluations of the real function, the prediction and its MSE
x = np.linspace(0,10,1000)

# Loop correlation models
corrs = (correxp1, correxp2, corrlin, corrcubic)
colors = ('b', 'g', 'y', 'm')
for k in range(len(corrs)):
    
    # Instanciate a Gaussian Process Model with the k-th correlation model
    if corrs[k] == corrlin or corrs[k] == corrcubic:
        aGaussianProcessModel = GaussianProcessModel(corr=corrs[k], theta0=1e-2, thetaL=1e-4, thetaU=1e-1, random_start=100)
    else:
        aGaussianProcessModel = GaussianProcessModel(corr=corrs[k], theta0=1e-2, thetaL=1e-4, thetaU=1e+1, random_start=100)
    
    # Fit to data using Maximum Likelihood Estimation of the parameters
    aGaussianProcessModel.fit(X, Y)
    
    # Make the prediction on the meshed x-axis (ask for MSE as well)
    y, MSE = aGaussianProcessModel.predict(x, eval_MSE=True)
    sigma = np.sqrt(MSE)
    
    # Compute the score function on a grid of the autocorrelation parameter space
    theta_values = np.logspace(np.log10(aGaussianProcessModel.thetaL[0,0]), np.log10(aGaussianProcessModel.thetaU[0,0]), 100)
    score_values = []
    for t in theta_values:
        score_values.append(aGaussianProcessModel.score(theta=t)[0])
    
    fig = pl.figure()
    
    # Plot the function, the prediction and the 95% confidence interval based on the MSE
    ax = fig.add_subplot(211)
    pl.plot(x, f(x), 'r:', label=u'$f(x) = x\,\sin(x)$')
    pl.plot(X, Y, 'r.', markersize=10, label=u'Observations')
    pl.plot(x, y, colors[k]+'-', label=u'Prediction (%s)' % corrs[k].__name__)
    pl.fill(np.concatenate([x, x[::-1]]), np.concatenate([y - 1.9600 * sigma, (y + 1.9600 * sigma)[::-1]]), alpha=.5, fc=colors[k], ec='None', label='95% confidence interval')
    pl.xlabel('$x$')
    pl.ylabel('$f(x)$')
    pl.ylim(-10, 20)
    pl.legend(loc='upper left')
    
    # Plot the score function
    ax = fig.add_subplot(212)
    pl.plot(theta_values, score_values, colors[k]+'-')
    pl.xlabel(u'$\\theta$')
    pl.ylabel(u'Score')
    pl.xscale('log')
    pl.xlim(aGaussianProcessModel.thetaL[0,0], aGaussianProcessModel.thetaU[0,0])

pl.show()
