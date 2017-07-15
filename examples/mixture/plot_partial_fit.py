
"""
===============================================
Bayesian Gaussian Mixture Model Fitting Methods
===============================================

Plot the confidence ellipsoids of a mixture of two Gaussians
obtained with Variational Inference (``BayesianGaussianMixture`` class models
with a Dirichlet process prior), using the ``fit()`` and ``partial_fit()``
methods.

The ``fit()`` method uses the entire dataset, either all at once or in
minibatches, while ``partial_fit()`` can be used for any number of points
equal to or greater than the number of componentsIn this case, we use
partial_fit each update on the minimum number of points necessary for fitting,
to simulate online learning
"""
# Authors: Joshua Engelman <j.aaron.engelman@gmail.com>
# License: BSD 3 clause

import itertools
import numpy as np
from scipy import linalg
import matplotlib.pyplot as plt
import matplotlib as mpl

from sklearn.mixture import BayesianGaussianMixture

print(__doc__)

color_iter = itertools.cycle(['navy', 'c', 'cornflowerblue', 'gold',
                              'darkorange'])


def plot_results(X, Y_, means, covariances, index, title):
    splot = plt.subplot(3, 1, 1 + index)
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
n_samples = 5000
n_components = 5

# Generate random sample, two components
np.random.seed(0)
X = np.zeros((n_samples, 2))
step = 4. * np.pi / n_samples

np.random.seed(0)
C = np.array([[0., -0.1], [1.7, .4]])
X = np.r_[np.dot(np.random.randn(n_samples, 2), C),
          .7 * np.random.randn(n_samples, 2) + np.array([-6, 3])]

plt.figure(figsize=(10, 10))
plt.subplots_adjust(bottom=.04, top=0.95, hspace=.2, wspace=.05,
                    left=.03, right=.97)

# Fit a Dirichlet process Gaussian mixture using the whole dataset
dpgmm = BayesianGaussianMixture(n_components=n_components,
                                random_state=1)

dpgmm.fit(X)
plot_results(X, dpgmm.predict(X), dpgmm.means_, dpgmm.covariances_, 0,
             'Fitting using full data')

# Fit a Dirichlet process Gaussian mixture using the whole dataset
dpgmm_minibatches = BayesianGaussianMixture(n_components=n_components,
                                            random_state=1,
                                            batch_size=100,
                                            n_jobs=-1)

dpgmm_minibatches.fit(X)
plot_results(X, dpgmm_minibatches.predict(X), dpgmm_minibatches.means_,
             dpgmm_minibatches.covariances_, 1,
             'Fitting using minibatches')

# Fit a Dirichlet process Gaussian mixture "online"

np.random.shuffle(X)
updates = np.array_split(X, np.floor(n_samples/n_components))
dpgmm_online = BayesianGaussianMixture(n_components=n_components,
                                       mean_precision_prior=1e-8,
                                       covariance_prior=1e1*np.eye(2),
                                       random_state=1)

scores = []
timestep = 0
for update in updates:
    timestep += 1
    dpgmm_online.partial_fit(update, lr=1./(timestep*1e-3))
    scores.append([timestep, dpgmm_online.score(X)])

plot_results(X, dpgmm_online.predict(X), dpgmm_online.means_,
             dpgmm_online.covariances_, 2,
             'Fitting online using partial_fit')

plt.show()
