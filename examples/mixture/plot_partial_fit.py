"""
===========================================
Bayesian Gaussian Mixture Model Ellipsoids
===========================================

Plot the confidence ellipsoids of a mixture of two Gaussians
obtained with Variational Inference (``BayesianGaussianMixture`` class models
with a Dirichlet process prior), using the ``fit()`` and ``partial_fit()``
methods.

The ``fit()`` method uses the entire dataset, while ``partial_fit()`` can be
used for either the full dataset or minibatches.  For streaming data or
datasets too large to fit into memory, ```partial_fit()``` gives a similar or
better fit much faster.
"""

import itertools

import numpy as np
from scipy import linalg
import matplotlib.pyplot as plt
import matplotlib as mpl

from sklearn import mixture

import time

color_iter = itertools.cycle(['navy', 'c', 'cornflowerblue', 'gold',
                              'darkorange'])


def plot_results(X, Y_, means, covariances, index, title):
    splot = plt.subplot(2, 1, 1 + index)
    for i, (mean, covar, color) in enumerate(zip(
            means, covariances, color_iter)):
        v, w = linalg.eigh(covar)
        v = 2. * np.sqrt(2.) * np.sqrt(v)
        u = w[0] / linalg.norm(w[0])
        # as the DP will not use every component it has access to
        # unless it needs it, we shouldn't plot the redundant
        # components.
        if not np.any(Y_ == i):
            continue
        plt.scatter(X[Y_ == i, 0], X[Y_ == i, 1], .8, color=color)

        # Plot an ellipse to show the Gaussian component
        angle = np.arctan(u[1] / u[0])
        angle = 180. * angle / np.pi  # convert to degrees
        ell = mpl.patches.Ellipse(mean, v[0], v[1], 180. + angle, color=color)
        ell.set_clip_box(splot.bbox)
        ell.set_alpha(0.5)
        splot.add_artist(ell)

    plt.xlim(-9., 5.)
    plt.ylim(-3., 6.)
    plt.xticks(())
    plt.yticks(())
    plt.title(title)


# Number of samples per component
n_samples = 10000

# Generate random sample, two components
np.random.seed(0)
C = np.array([[0., -0.1], [1.7, .4]])
X = np.r_[np.dot(np.random.randn(n_samples, 2), C),
          .7 * np.random.randn(n_samples, 2) + np.array([-6, 3])]

# Fit a Dirichlet process Gaussian mixture using five components
dpgmm = mixture.BayesianGaussianMixture(n_components=5,
                                        covariance_type='full',
                                        max_iter=200)
start = time.time()
dpgmm.fit(X)
full_runtime = time.time() - start

plot_results(X, dpgmm.predict(X), dpgmm.means_, dpgmm.covariances_, 0,
             'Fitting using full data')

dpgmm2 = mixture.BayesianGaussianMixture(n_components=5,
                                         covariance_type='full',
                                         max_iter=200, warm_start=True)
n_batches = 10
minibatches = np.array_split(X, n_batches)

start2 = time.time()
for batch in minibatches:
    dpgmm2.partial_fit(batch)

partial_runtime = (time.time() - start2)

plot_results(X, dpgmm2.predict(X), dpgmm2.means_, dpgmm2.covariances_, 1,
             'Fitting using batches')

print "fit() runtime: ", full_runtime
print "partial_fit() runtime: ", partial_runtime
plt.show()
