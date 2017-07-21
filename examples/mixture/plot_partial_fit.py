
"""
=================================================
Online Fitting of Bayesian Gaussian Mixture Model
=================================================

Plot the confidence ellipsoids of a mixture of two Gaussians
obtained with Variational Inference (``BayesianGaussianMixture`` class models
with a Dirichlet process prior), using the ``fit()`` and ``partial_fit()``
methods.

The ``fit()`` method uses the entire dataset, either all at once or in
minibatches, while ``partial_fit()`` can be used for any number of points
equal to or greater than the number of componentsI.  In this case, we use
partial_fit on batches of size 500, and restrict each batch to 10 iterations,
We see that we get nearly identical results in a fraction of the time.
The effect is more pronounced if we start with a larger number of components.
"""
# Authors: Joshua Engelman <j.aaron.engelman@gmail.com>
# License: BSD 3 clause

from itertools import cycle
from time import time

import numpy as np
from scipy import linalg

import matplotlib.pyplot as plt
import matplotlib as mpl

from sklearn.datasets import make_blobs
from sklearn.mixture import BayesianGaussianMixture
from sklearn.externals.six.moves import xrange

print(__doc__)


color_iter = cycle(['navy', 'c', 'cornflowerblue', 'gold',
                    'darkorange', 'brown', 'k', 'g', 'r', 'm'])


def plot_results(X, estimator, title, subplot):

    means = estimator.means_
    covariances = estimator.covariances_
    Y_ = estimator.predict(X)
    loss = estimator.score(X)

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
        subplot.scatter(X[Y_ == i, 0], X[Y_ == i, 1], .8, color=color)

        # Plot an ellipse to show the Gaussian component
        angle = np.arctan(u[1] / u[0])
        angle = 180. * angle / np.pi  # convert to degrees
        ell = mpl.patches.Ellipse(mean, v[0], v[1], 180. + angle, color=color)
        ell.set_clip_box(subplot.bbox)
        ell.set_alpha(0.5)
        subplot.add_artist(ell)

    subplot.set_xticks(())
    subplot.set_xlabel('Log likelihood: %f' % loss)
    subplot.set_yticks(())
    subplot.set_title(title)


# Number of samples per component
n_samples = 10000
batch_size = 500

# Generate random sample, two components
np.random.seed(0)

point = 5 * np.sqrt(2) / 2.
centers = [[-5, 0], [-point, point], [5, 0], [point, point], [0, 0],
           [5, 0], [point, -point], [0, -5], [-point, -point], [0, 5]]
X, labels = make_blobs(n_samples=n_samples, centers=centers,
                       cluster_std=1.1)

updates = np.array_split(X, np.ceil(n_samples / (batch_size)))


f, axarr = plt.subplots(2, 2, figsize=(10, 10))

n_components = [20, 12]

for i in xrange(len(n_components)):
    # Fit a Dirichlet process Gaussian Mixture using the full dataset
    dpgmm = BayesianGaussianMixture(n_components=n_components[i],
                                    random_state=1,
                                    max_iter=500)

    start = time()
    dpgmm.fit(X)
    runtime = time() - start

    title = '''fit() with %i components,
                runtime %.3f s''' % (n_components[i], runtime)
    plot_results(X, dpgmm, title, axarr[0, i])

    # Fit a Dirichlet process Gaussian mixture using minibatches
    dpgmm_online = BayesianGaussianMixture(n_components=n_components[i],
                                           random_state=1,
                                           max_iter=10)

    start2 = time()
    for update in updates:
        dpgmm_online.partial_fit(update)

    runtime2 = time() - start2

    title = 'partial_fit() with 20 components, runtime %.3f s' % runtime2
    plot_results(X, dpgmm_online, title, axarr[1, i])


plt.show()
