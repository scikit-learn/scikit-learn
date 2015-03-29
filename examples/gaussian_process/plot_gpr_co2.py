"""Gaussian process regression (GPR) on Mauna Loa CO2 data. """
print __doc__

# Authors: Jan Hendrik Metzen <jhm@informatik.uni-bremen.de>
#
# License: BSD 3 clause

import numpy as np
import statsmodels.api as sm   # XXX: Upload data on mldata

from matplotlib import pyplot as plt

from sklearn.gaussian_process import GaussianProcessRegressor
from sklearn.gaussian_process.kernels \
    import RBF, Kernel, WhiteKernel, RationalQuadratic, ExpSineSquared

data = sm.datasets.get_rdataset("co2").data
X = np.array(data.time)[:, np.newaxis]
y = np.array(data.co2)
y_mean = y.mean()

# Kernel with parameters given in GPML book
k1 = 66.0**2 * RBF(l=67.0)  # long term smooth rising trend
k2 = 2.4**2 * RBF(l=90.0) * ExpSineSquared(l=1.3, p=1.0)  # seasonal component
k3 = 0.66**2 * RationalQuadratic(alpha=0.78, l=1.2)  # medium term irregularit.
k4 = 0.18**2 * RBF(l=0.134) + WhiteKernel(c=0.19**2) # noise terms
kernel_gpml = k1 + k2 + k3 + k4

gp = GaussianProcessRegressor(kernel=kernel_gpml, y_err=0, optimizer=None)
gp.fit(X, y - y_mean)

print "GPML kernel: %s" % gp.kernel_
print "Log-marginal-likelihood: %.3f" % gp.log_marginal_likelihood(gp.theta_)

# Kernel with optimized parameters
k1 = 50.0**2 * RBF(l=50.0) # long term smooth rising trend
k2 = 2.0**2 * RBF(l=100.0) \
	* ExpSineSquared(l=1.0, p=1.0, p_bounds="fixed") # seasonal component
k3 = 0.5**2 * RationalQuadratic(alpha=1.0, l=1.0) # medium term irregularities
k4 = 0.1**2 * RBF(l=0.1) + WhiteKernel(c=0.1**2,
									   c_bounds=(1e-3, np.inf))  # noise terms
kernel = k1 + k2 + k3 + k4

gp = GaussianProcessRegressor(kernel=kernel, y_err=0)
gp.fit(X, y - y_mean)

print "\nLearned kernel: %s" % gp.kernel_
print "Log-marginal-likelihood: %.3f" % gp.log_marginal_likelihood(gp.theta_)

X_ = np.linspace(X.min(), X.max() + 30, 1000)[:, np.newaxis]
y_pred, y_std = gp.predict(X_, return_std=True)
y_pred += y_mean

# Illustration
plt.scatter(X, y, c='k')
plt.plot(X_, y_pred)
plt.fill_between(X_[:, 0], y_pred - y_std, y_pred + y_std,
                 alpha=0.5, color='k')
plt.xlim(X_.min(), X_.max())
plt.xlabel("Year")
plt.ylabel(r"CO$_2$ in ppm")
plt.title(r"Atmospheric CO$_2$ concentration at Mauna Loa")
plt.tight_layout()
plt.show()
