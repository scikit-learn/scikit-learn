#!/usr/bin/python
# -*- coding: utf-8 -*-

r"""
=========================================================
Gaussian Processes regression: non-stationary correlation
=========================================================

A simple one-dimensional regression task which is complicated by a
discontinuity at x=0. Because of this discontinuity, the global GP
hyperparameters are chosen such that the correlation drops quickly with
increasing distance of the datapoints. This is inappropriate for all pairs of
datapoints which are not separated by the discontinuity.

A non-stationary correlation model, which assigns shorter length-scales to
regions close to the discontinuity can reduce this issue. In this script,
a local length-scale model is provided and only a global modifier of the
length scales is learned. This non-stationary model generalizes more globally
and obtains a considerably higher posterior probability.
"""
print(__doc__)

# Author: Jan Hendrik Metzen <jhm@informatik.uni-bremen.de>
# Licence: BSD 3 clause

import numpy as np
from scipy.special import gamma, kv
from sklearn.gaussian_process import GaussianProcess
from sklearn.gaussian_process.correlation_models \
    import MACHINE_EPSILON, l1_cross_differences
from matplotlib import pyplot as pl

np.random.seed(0)


def f(x):
    """The function to predict.

    The function has a discontinuity at x=0.
    """
    return np.where(x[:, 0] > 0, np.sin(x[:, 0]), np.cos(-x[:, 0])) - 2


#----------------------------------------------------------------------
# Correlation model

def length_scale(X):
    """ Predetermined length-scales for the test problem.

    Length-scales are smaller close to the discontinuity at x=0
    """
    return 10**np.minimum(0, 5 * X[:, 0] ** 2 - 1.0)


class NonStationaryCorrelation(object):
    """ Non-stationary correlation model with predetermined length-scales."""

    def fit(self, X, nugget=10. * MACHINE_EPSILON):
        self.X = X
        self.nugget = nugget
        self.n_samples = self.X.shape[0]

        # Calculate array with shape (n_eval, n_features) giving the
        # componentwise distances between locations x and x' at which the
        # correlation model should be evaluated.
        self.D, self.ij = l1_cross_differences(self.X)

        # Calculate length scales
        self.l_train = length_scale(self.X)

    def __call__(self, theta, X=None):
        # Prepare distances and length scale information for any pair of
        # datapoints, whose correlation shall be computed
        if X is not None:
            # Get pairwise componentwise L1-differences to the input training
            # set
            d = X[:, np.newaxis, :] - self.X[np.newaxis, :, :]
            d = d.reshape((-1, X.shape[1]))
            # Calculate length scales
            l_query = length_scale(X)
            l = np.transpose([np.tile(self.l_train, len(l_query)),
                              np.repeat(l_query, len(self.l_train))])
        else:
            # No external datapoints given; auto-correlation of training set
            # is used instead
            d = self.D
            l = self.l_train[self.ij]

        # Compute general Matern kernel for nu=1.5
        nu = 1.5
        if d.ndim > 1 and theta.size == d.ndim:
            activation = np.sum(theta.reshape(1, d.ndim) * d ** 2, axis=1)
        else:
            activation = theta[0] * np.sum(d ** 2, axis=1)
        tmp = 0.5*(l**2).sum(1)
        tmp2 = 2*np.sqrt(nu * activation / tmp)
        r = np.sqrt(l[:, 0]) * np.sqrt(l[:, 1]) / (gamma(nu) * 2**(nu - 1))
        r /= np.sqrt(tmp)
        r *= tmp2**nu * kv(nu, tmp2)

        # Convert correlations to 2d matrix
        if X is not None:
            return r.reshape(-1, self.n_samples)
        else:  # exploit symmetry of auto-correlation
            R = np.eye(self.n_samples) * (1. + self.nugget)
            R[self.ij[:, 0], self.ij[:, 1]] = r
            R[self.ij[:, 1], self.ij[:, 0]] = r
            return R

    def log_prior(self, theta):
        # Just use a flat prior
        return 0

#----------------------------------------------------------------------
# Actual test data
X = np.random.random(50)[:, None] * 4 - 2

# Observations
y = f(X).ravel()

# Mesh the input space for evaluations of the real function, the prediction and
# its MSE
x = np.atleast_2d(np.linspace(-2, 2, 1000)).T

# Instanciate one Gaussian Process model for the stationary Matern kernel and
# one for the non-stationary one
gp_stationary = \
    GaussianProcess(corr='matern_1.5', theta0=1e0, thetaL=1e-2, thetaU=1e+2,
                    random_start=100)
gp_non_stationary = \
    GaussianProcess(corr=NonStationaryCorrelation(),
                    theta0=1e0, thetaL=1e-2, thetaU=1e+2,
                    random_start=100)

# Fit to data using Maximum Likelihood Estimation of the parameters
gp_stationary.fit(X, y)
gp_non_stationary.fit(X, y)
print("Theta:\n\tStationary: {:.3f} \t Non-stationary: {:.3f}"
      .format(gp_stationary.theta_[0], gp_non_stationary.theta_[0]))
print("Posterior probability (negative, average, log):\n\t"
      "Stationary: {:.5f} \t Non-stationary: {:.5f}"
      .format(gp_stationary.posterior_function_value_,
              gp_non_stationary.posterior_function_value_))

# Plot predictions
for title, gp in [("stationary", gp_stationary),
                  ("non-stationary", gp_non_stationary)]:
    # Make the prediction on the meshed x-axis (ask for MSE as well)
    y_pred, MSE = gp.predict(x, eval_MSE=True)
    sigma = np.sqrt(MSE)

    # Plot the function, the prediction and the 95% confidence interval
    # based on the MSE
    fig = pl.figure()
    pl.plot(x, f(x), 'r:',
            label=u'$f(x) = -2 + sin(x)$ if $x > 0$ else $-2 + cos(-x)$')
    pl.plot(X, y, 'r.', markersize=10, label=u'Observations')
    pl.plot(x, y_pred, 'b-', label=u'Prediction')
    pl.fill(np.concatenate([x, x[::-1]]),
            np.concatenate([y_pred - 1.9600 * sigma,
                           (y_pred + 1.9600 * sigma)[::-1]]),
            alpha=.5, fc='b', ec='None', label='95% confidence interval')
    if gp == gp_non_stationary:
        pl.plot(x, length_scale(x), label="Length scales")
    pl.xlabel('$x$')
    pl.ylabel('$f(x)$')
    pl.ylim(-5, 5)
    pl.legend(loc='upper left')
    pl.title(title)

pl.show()
