
import numpy as np
import pylab

from sklearn.gaussian_process import GaussianProcessRegression
from sklearn.gaussian_process.kernels import RBF

np.random.seed(0)


# Specify Gaussian Process
kernel = (1e-10, 1.0, 100) * RBF(param_space=(1e-10, 1.0, None))
gp = GaussianProcessRegression(kernel=kernel)

# Plot prior
pylab.figure(0, figsize=(8, 8))
pylab.subplot(2, 1, 1)
X_ = np.linspace(0, 5, 100)
y_mean, y_cov = gp.predict(X_[:, None], return_cov=True)
pylab.plot(X_, y_mean, 'k', lw=3, zorder=9)
pylab.fill_between(X_, y_mean - np.sqrt(np.diag(y_cov)),
                   y_mean + np.sqrt(np.diag(y_cov)),
                   alpha=0.5, color='k')
y_samples = gp.sample(X_[:, None], 10)
pylab.plot(X_, y_samples, color='b', lw=1)
pylab.xlim(0, 5)
pylab.ylim(-3, 3)
pylab.title("Prior")

# Generate data and fit GP
X = np.random.uniform(0, 5, 10)[:, None]
y = np.sin((X[:, 0] - 2.5) ** 2)
gp.fit(X, y)

# Plot posterior
pylab.subplot(2, 1, 2)
X_ = np.linspace(0, 5, 100)
y_mean, y_cov = gp.predict(X_[:, None], return_cov=True)
pylab.plot(X_, y_mean, 'k', lw=3, zorder=9)
pylab.fill_between(X_, y_mean - np.sqrt(np.diag(y_cov)),
                   y_mean + np.sqrt(np.diag(y_cov)),
                   alpha=0.5, color='k')

y_samples = gp.sample(X_[:, None], 10)
pylab.plot(X_, y_samples, color='b', lw=1)
pylab.scatter(X[:, 0], y, c='r', s=50, zorder=10)
pylab.xlim(0, 5)
pylab.ylim(-3, 3)
pylab.title("Posterior")
pylab.tight_layout()
pylab.show()
