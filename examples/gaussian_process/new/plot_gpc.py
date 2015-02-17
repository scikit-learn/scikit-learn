
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
pylab.plot(X_, gp_fix.predict_proba(X_[:, None]), 'r', label="Fixed")
pylab.plot(X_, gp_opt.predict_proba(X_[:, None]), 'b', label="Optimized")
pylab.legend(loc="best")

# Plot LML landscape
pylab.figure(1)
theta0 = np.logspace(0, 8, 50)
theta1 = np.logspace(-1, 1, 49)
Theta0, Theta1 = np.meshgrid(theta0, theta1)
LML = [[gp_opt.log_marginal_likelihood([Theta0[i, j], Theta1[i, j]])
        for i in range(Theta0.shape[0])] for j in range(Theta0.shape[1])]
LML = np.array(LML).T
pylab.pcolor(Theta0, Theta1, LML)
pylab.xscale("log")
pylab.yscale("log")
pylab.colorbar()
pylab.show()
