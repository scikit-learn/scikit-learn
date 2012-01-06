#!/usr/bin/python
# -*- coding: utf-8 -*-

"""
========================================================================
Gaussian Processes regression: goodness-of-fit on the 'diabetes' dataset
========================================================================

This example consists in fitting a Gaussian Process model onto the diabetes
dataset.

The correlation parameters are determined by means of maximum likelihood
estimation (MLE). An anisotropic squared exponential correlation model with a
constant regression model are assumed. We also used a nugget = 1e-2 in order to
account for the (strong) noise in the targets.

We compute then compute a cross-validation estimate of the coefficient of
determination (R2) without reperforming MLE, using the set of correlation
parameters found on the whole dataset.
"""
print __doc__

# Author: Vincent Dubourg <vincent.dubourg@gmail.com>
# License: BSD style

from sklearn import datasets
from sklearn.gaussian_process import GaussianProcess
from sklearn.cross_validation import cross_val_score, KFold

# Load the dataset from scikit's data sets
diabetes = datasets.load_diabetes()
X, y = diabetes.data, diabetes.target

# Instanciate a GP model
gp = GaussianProcess(regr='constant', corr='absolute_exponential',
                     theta0=[1e-4] * 10, thetaL=[1e-12] * 10,
                     thetaU=[1e-2] * 10, nugget=1e-2, optimizer='Welch')

# Fit the GP model to the data performing maximum likelihood estimation
gp.fit(X, y)

# Deactivate maximum likelihood estimation for the cross-validation loop
gp.theta0 = gp.theta  # Given correlation parameter = MLE
gp.thetaL, gp.thetaU = None, None  # None bounds deactivate MLE

# Perform a cross-validation estimate of the coefficient of determination using
# the cross_validation module using all CPUs available on the machine
K = 20  # folds
R2 = cross_val_score(gp, X, y=y, cv=KFold(y.size, K), n_jobs=1).mean()
print("The %d-Folds estimate of the coefficient of determination is R2 = %s"
    % (K, R2))
