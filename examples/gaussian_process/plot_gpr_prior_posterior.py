"""Gaussian process regression (GPR) prior and posterior

This example illustrates the prior and posterior of a GPR with different
kernels. Mean, standard deviation, and 10 samples are shown for both prior
and posterior.
"""
print __doc__

# Authors: Jan Hendrik Metzen <jhm@informatik.uni-bremen.de>
#
# License: BSD 3 clause

import numpy as np

from matplotlib import pyplot as plt

from sklearn.gaussian_process import GaussianProcessRegressor
from sklearn.gaussian_process.kernels \
    import RBF, RationalQuadratic, ExpSineSquared, DotProduct


kernels = [1.0 * RBF(1.0, 1e-1, 10.0),
           1.0 * RationalQuadratic((0.1, 1.0)),
           1.0 * ExpSineSquared((1.0, 3.0), (0.1, 1.0), (10.0, 10.0)),
           (0.01, 0.1, 10.0) * (DotProduct(1.0, 0.0, 10.0) ** 2)]

for fig_index, kernel in enumerate(kernels):
    if fig_index > 3: continue

    # Specify Gaussian Process
    gp = GaussianProcessRegressor(kernel=kernel)

    # Plot prior
    plt.figure(fig_index, figsize=(8, 8))
    plt.subplot(2, 1, 1)
    X_ = np.linspace(0, 5, 100)
    y_mean, y_cov = gp.predict(X_[:, np.newaxis], return_cov=True)
    plt.plot(X_, y_mean, 'k', lw=3, zorder=9)
    plt.fill_between(X_, y_mean - np.sqrt(np.diag(y_cov)),
                     y_mean + np.sqrt(np.diag(y_cov)),
                     alpha=0.5, color='k')
    y_samples = gp.sample_y(X_[:, np.newaxis], 10)
    plt.plot(X_, y_samples, color='b', lw=1)
    plt.xlim(0, 5)
    plt.ylim(-3, 3)
    plt.title("Prior (kernel:  %s)" % kernel)

    # Generate data and fit GP
    rng = np.random.RandomState(4)
    X = rng.uniform(0, 5, 10)[:, np.newaxis]
    y = np.sin((X[:, 0] - 2.5) ** 2)
    gp.fit(X, y)

    # Plot posterior
    plt.subplot(2, 1, 2)
    X_ = np.linspace(0, 5, 100)
    y_mean, y_cov = gp.predict(X_[:, np.newaxis], return_cov=True)
    plt.plot(X_, y_mean, 'k', lw=3, zorder=9)
    plt.fill_between(X_, y_mean - np.sqrt(np.diag(y_cov)),
                     y_mean + np.sqrt(np.diag(y_cov)),
                     alpha=0.5, color='k')

    y_samples = gp.sample_y(X_[:, np.newaxis], 10)
    plt.plot(X_, y_samples, color='b', lw=1)
    plt.scatter(X[:, 0], y, c='r', s=50, zorder=10)
    plt.xlim(0, 5)
    plt.ylim(-3, 3)
    plt.title("Posterior (kernel: %s)" % gp.kernel_)
    plt.tight_layout()

plt.show()
