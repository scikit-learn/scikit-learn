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

from matplotlib import pyplot as plt

from sklearn.gaussian_process import GaussianProcessClassifier
from sklearn.gaussian_process.kernels import RBF


# Generate data
rng = np.random.RandomState(0)
X = rng.uniform(0, 5, 50)[:, np.newaxis]
y = np.array(np.sin((X[:, 0] - 2.5) ** 2) > 0.0, dtype=int)

# Specify Gaussian Processes with fixed and optimized hyperparameters
kernel_fix = 1.0 * RBF(param_space=1.0)
gp_fix = GaussianProcessClassifier(kernel=kernel_fix, optimizer=None).fit(X, y)

kernel_opt = 1.0 * RBF(1.0)
gp_opt = GaussianProcessClassifier(kernel=kernel_opt).fit(X, y)

print "Log Marginal Likelihood (initial): %.3f" % \
    gp_fix.log_marginal_likelihood(gp_fix.theta_)
print "Log Marginal Likelihood (optimized): %.3f" % \
    gp_fix.log_marginal_likelihood(gp_opt.theta_)


# Plot posteriors
plt.figure(0)
plt.scatter(X[:, 0], y)
X_ = np.linspace(0, 5, 100)
plt.plot(X_, gp_fix.predict_proba(X_[:, np.newaxis]), 'r',
         label="Initial kernel: %s" % kernel_fix)
plt.plot(X_, gp_opt.predict_proba(X_[:, np.newaxis]), 'b',
         label="Optimized kernel: %s" % kernel_opt)
plt.legend(loc="best")
plt.xlabel("Feature")
plt.ylabel("Class")

# Plot LML landscape
plt.figure(1)
theta0 = np.logspace(0, 8, 30)
theta1 = np.logspace(-1, 1, 29)
Theta0, Theta1 = np.meshgrid(theta0, theta1)
LML = [[gp_opt.log_marginal_likelihood([Theta0[i, j], Theta1[i, j]])
        for i in range(Theta0.shape[0])] for j in range(Theta0.shape[1])]
LML = np.array(LML).T
plt.pcolor(Theta0, Theta1, LML)
plt.xscale("log")
plt.yscale("log")
plt.colorbar()
plt.xlabel("Magnitude")
plt.ylabel("Length-scale")
plt.title("Log-marginal-likelihood")

plt.show()
