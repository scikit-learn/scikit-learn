"""
==================
Gaussian mixture clustering
==================

Demonstration of Gaussian mixture models for clustering.

See :ref:`gmm` for more information on the estimator.

Plots the fitted Gaussian mixture distributions using a
variety of GMM classifiers on the toy data.

Compares Gaussian mixture with spherical, diagonal, full, and tied covariance.
"""
print(__doc__)

# Author: Wei Xue <xuewei4d@gmail.com>
# License: BSD 3 clause

# $Id$

import matplotlib.pyplot as plt
import matplotlib as mpl
import numpy as np
from sklearn.datasets.samples_generator import make_spd_matrix
from sklearn.mixture import GaussianMixture as GM


rng = np.random.RandomState(42)
n_samples = 100
n_features = 2
n_components = 3


def make_ellipses(gm, ax):
    for n, color in enumerate('rgb'):
        if gm.covariance_type == 'full':
            covars = gm.covars_[n]
        elif gm.covariance_type == 'tied':
            covars = gm.covars_
        elif gm.covariance_type == 'diag':
            covars = np.diag(gm.covars_[n])
        elif gm.covariance_type == 'spherical':
            covars = np.eye(n_features) * gm.covars_[n]
        v, w = np.linalg.eigh(covars)
        u = w[0] / np.linalg.norm(w[0])
        angle = np.arctan2(u[1], u[0])
        angle = 180 * angle / np.pi  # convert to degrees
        v *= 2
        ell = mpl.patches.Ellipse(gm.means_[n], v[0], v[1],
                                  180 + angle, color=color)
        ell.set_clip_box(ax.bbox)
        ell.set_alpha(0.5)
        ax.add_artist(ell)

means = [rng.uniform(-1, 1, n_features) * 5 for _ in range(n_components)]
covars = [make_spd_matrix(n_features, rng) for _ in range(n_components)]
weights = rng.rand(n_components)
weights = weights / np.sum(weights)

n_samples_components = np.round(weights * n_samples).astype(int)
X = np.vstack([rng.multivariate_normal(
    means[j], covars[j], np.round(weights[j] * n_samples).astype(int))
    for j in range(n_components)])

gms = dict((covar_type, GM(
    n_components=n_components, covariance_type=covar_type,
    init_params='kmeans', n_iter=100, reg_covar=0))
    for covar_type in ['spherical', 'diag', 'tied', 'full'])

plt.figure()
for index, (name, gm) in enumerate(gms.items()):
    gm.fit(X)
    h = plt.subplot(2, len(gms) / 2, index + 1)
    h.scatter(X[:, 0], X[:, 1], s=5)
    make_ellipses(gm, h)
    plt.title(name)

plt.legend(loc='lower right', prop=dict(size=12))
plt.show()
