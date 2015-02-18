"""Gaussian process classification (GPC)

This example illustrates both prediction of the prior GPC and the posterior
GPC. While the posterior model has a considerably larger
log-marginal-likelihood, the generated predictions are not optimal. This
is caused by the Laplace approximations used internally by GPC.
"""
print __doc__

# Authors: Jan Hendrik Metzen <jhm@informatik.uni-bremen.de>
#
# License: BSD 3 clause

import numpy as np

from sklearn.gaussian_process import GaussianProcessClassification
from sklearn.gaussian_process.kernels import RBF

np.random.seed(0)


# Generate data
X = np.random.uniform(0, 5, 50)[:, None]
y = np.array(np.sin((X[:, 0] - 2.5) ** 2) > 0.0, dtype=int)

# Specify Gaussian Processes with fixed and optimized hyperparameters
kernel_fix = 4.0 * RBF(param_space=[1.0])
gp_fix = GaussianProcessClassification(kernel=kernel_fix).fit(X, y)

kernel_opt = (1e-10, 1.0, 100) * RBF(param_space=(1e-10, 1, 10))
gp_opt = GaussianProcessClassification(kernel=kernel_opt).fit(X, y)

print "Log Marginal Likelihood (initial): %.3f" % \
    gp_fix.log_marginal_likelihood(gp_fix.theta_)
print "Log Marginal Likelihood (optimized): %.3f" % \
    gp_fix.log_marginal_likelihood(gp_opt.theta_)


# Plot posteriors
import pylab
pylab.figure(0)
pylab.scatter(X[:, 0], y)
X_ = np.linspace(0, 5, 100)
pylab.plot(X_, gp_fix.predict_proba(X_[:, None]), 'r',
           label="Initial kernel: %s" % kernel_fix)
pylab.plot(X_, gp_opt.predict_proba(X_[:, None]), 'b',
           label="Optimized kernel: %s" % kernel_opt)
pylab.legend(loc="best")
pylab.xlabel("Feature")
pylab.ylabel("Class")

# Plot LML landscape
pylab.figure(1)
theta0 = np.logspace(0, 8, 30)
theta1 = np.logspace(-1, 1, 29)
Theta0, Theta1 = np.meshgrid(theta0, theta1)
LML = [[gp_opt.log_marginal_likelihood([Theta0[i, j], Theta1[i, j]])
        for i in range(Theta0.shape[0])] for j in range(Theta0.shape[1])]
LML = np.array(LML).T
pylab.pcolor(Theta0, Theta1, LML)
pylab.xscale("log")
pylab.yscale("log")
pylab.colorbar()
pylab.xlabel("Magnitude")
pylab.ylabel("Length-scale")
pylab.title("Log-marginal-likelihood")

pylab.show()
