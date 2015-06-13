"""
==================================================================================
Gaussian Mixture, Bayesian Gaussian Mixture and Dirichlet Process Gaussian
Mixture
==================================================================================

Plot the confidence ellipsoids of a mixture of two Gaussians with
GaussianMixture, Bayesia Gaussian Mixture and Dirichlet Process Gaussian
Mixture

The three models have access to five components with which to fit the
data. Note that the GM model will necessarily use all five components
while BGM and DPGM model will effectively only use as many as are needed for
a good fit. This is a property of the Dirichlet distribution prior and
Dirichlet process prior. Here we can see that the GM model splits some
components arbitrarily, because it is trying to fit too many components, while
BGM and DPGM adapt it number of state automatically.

This example doesn't show it, as we're in a low-dimensional space, but
another advantage of the Dirichlet process model is that it can fit
full covariance matrices effectively even when there are less examples
per cluster than there are dimensions in the data, due to
regularization properties of the inference algorithm.
"""
import itertools

import numpy as np
from scipy import linalg
import matplotlib.pyplot as plt
import matplotlib as mpl

from sklearn import mixture
from sklearn.datasets.samples_generator import make_spd_matrix


# Number of samples per component
n_samples = 500

# Generate random sample, three components
rng = np.random.RandomState(42)
n_features = 2
n_components = 3

means = [rng.uniform(-1, 1, n_features) * 5 for _ in range(n_components)]
covars = [make_spd_matrix(n_features, rng) for _ in range(n_components)]
weights = rng.rand(n_components)
weights = weights / np.sum(weights)

n_samples_components = np.round(weights * n_samples).astype(int)
X = np.vstack([rng.multivariate_normal(
    means[j], covars[j], np.round(weights[j] * n_samples).astype(int))
    for j in range(n_components)])

# Fit a mixture of Gaussians with EM
gm = mixture.GaussianMixture(n_components=5, covariance_type='full',
                             reg_covar=0, init_params='kmeans')
gm.fit(X)

# Fit a mixture of Bayesian Gaussians with variational inference
bgm = mixture.BayesianGaussianMixture(n_components=5, precision_type='full',
                                      reg_covar=0, init_params='random',
                                      alpha_prior=1e-2)
bgm.fit(X)

# Fit a Dirichlet process mixture of Gaussians with variational inference
dpgm = mixture.DirichletProcessGaussianMixture(n_components=5,
                                               precision_type='full',
                                               init_params='random',
                                               gamma_prior=1e-2)
dpgm.fit(X)

color_iter = itertools.cycle(['r', 'g', 'b', 'c', 'm'])

for i, (clf, title) in enumerate([(gm, 'GM'),
                                  (bgm, 'BGM'),
                                  (dpgm, 'DPGM')]):
    splot = plt.subplot(3, 1, 1 + i)
    Y_ = clf.predict(X)
    for i, (mean, covar, color) in enumerate(zip(
            clf.means_, clf.covars_, color_iter)):
        v, w = linalg.eigh(covar)
        u = w[0] / linalg.norm(w[0])
        # as the DP will not use every component it has access to
        # unless it needs it, we shouldn't plot the redundant
        # components.
        if not np.any(Y_ == i):
            continue
        plt.scatter(X[Y_ == i, 0], X[Y_ == i, 1], .8, color=color)

        # Plot an ellipse to show the Gaussian component
        angle = np.arctan(u[1] / u[0])
        angle = 180 * angle / np.pi  # convert to degrees
        ell = mpl.patches.Ellipse(mean, v[0], v[1], 180 + angle, color=color)
        ell.set_clip_box(splot.bbox)
        ell.set_alpha(0.5)
        splot.add_artist(ell)

    plt.xticks(())
    plt.yticks(())
    plt.title(title)

plt.show()
