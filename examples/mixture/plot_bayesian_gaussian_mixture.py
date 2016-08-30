"""
======================================================
Bayesian Gaussian Mixture Concentration Prior Analysis
======================================================

Plot the resulting ellipsoids of a mixture of three Gaussians with
variational Bayesian Gaussian Mixture for three different values on the
prior the dirichlet concentration.

For all models, the Variationnal Bayesian Gaussian Mixture adapts its number of
mixture automatically. The parameter `dirichlet_concentration_prior` has a
direct link with the resulting number of components. Specifying a high value of
`dirichlet_concentration_prior` leads more often to uniformly-sized mixture
components, while specifying small (under 0.1) values will lead to some mixture
components getting almost all the points while most mixture components will be
centered on just a few of the remaining points.
"""
# Author: Thierry Guillemot <thierry.guillemot.work@gmail.com>
# License: BSD 3 clause

import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec

from sklearn.mixture import BayesianGaussianMixture

print(__doc__)


def plot_ellipses(ax, weights, means, covars):
    for n in range(means.shape[0]):
        v, w = np.linalg.eigh(covars[n][:2, :2])
        u = w[0] / np.linalg.norm(w[0])
        angle = np.arctan2(u[1], u[0])
        angle = 180 * angle / np.pi  # convert to degrees
        v = 2 * np.sqrt(2) * np.sqrt(v)
        ell = mpl.patches.Ellipse(means[n, :2], v[0], v[1], 180 + angle)
        ell.set_clip_box(ax.bbox)
        ell.set_alpha(weights[n])
        ax.add_artist(ell)


def plot_results(ax1, ax2, estimator, dirichlet_concentration_prior, X, y, plot_title=False):
    estimator.dirichlet_concentration_prior = dirichlet_concentration_prior
    estimator.fit(X)
    ax1.set_title("Bayesian Gaussian Mixture for "
                  r"$dc_0=%.1e$" % dirichlet_concentration_prior)
    # ax1.axis('equal')
    ax1.scatter(X[:, 0], X[:, 1], s=5, marker='o', color=colors[y], alpha=0.8)
    ax1.set_xlim(-2., 2.)
    ax1.set_ylim(-3., 3.)
    ax1.set_xticks(())
    ax1.set_yticks(())
    plot_ellipses(ax1, estimator.weights_, estimator.means_,
                  estimator.covariances_)

    ax2.get_xaxis().set_tick_params(direction='out')
    ax2.yaxis.grid(True, alpha=0.7)
    for k, w in enumerate(estimator.weights_):
        ax2.bar(k - .45, w, width=0.9, color='royalblue', zorder=3)
        ax2.text(k, w + 0.007, "%.1f%%" % (w * 100.),
                 horizontalalignment='center')
    ax2.set_xlim(-.6, 2 * n_components - .4)
    ax2.set_ylim(0., 1.1)
    ax2.tick_params(axis='y', which='both', left='off',
                    right='off', labelleft='off')
    ax2.tick_params(axis='x', which='both', top='off')

    if plot_title:
        ax1.set_ylabel('Estimated Mixtures')
        ax2.set_ylabel('Weight of each component')

# Parameters
random_state = 2
n_components, n_features = 3, 2
colors = np.array(['mediumseagreen', 'royalblue', 'r', 'gold',
                   'orchid', 'indigo', 'darkcyan', 'tomato'])
dirichlet_concentration_prior = np.logspace(-3, 3, 3)
covars = np.array([[[.7, .0], [.0, .1]],
                   [[.5, .0], [.0, .1]],
                   [[.5, .0], [.0, .1]]])
samples = np.array([200, 500, 200])
means = np.array([[.0, -.70],
                  [.0, .0],
                  [.0, .70]])


# Here we put beta_prior to 0.8 to minimize the influence of the prior for this
# dataset
estimator = BayesianGaussianMixture(n_components=2 * n_components,
                                    init_params='random', max_iter=1500,
                                    mean_precision_prior=.8, tol=1e-9,
                                    random_state=random_state)

# Generate data
rng = np.random.RandomState(random_state)
X = np.vstack([
    rng.multivariate_normal(means[j], covars[j], samples[j])
    for j in range(n_components)])
y = np.concatenate([j * np.ones(samples[j], dtype=int)
                    for j in range(n_components)])

# Plot Results
plt.figure(figsize=(4.7 * 3, 8))
plt.subplots_adjust(bottom=.04, top=0.95, hspace=.05, wspace=.05,
                    left=.03, right=.97)

gs = gridspec.GridSpec(3, len(dirichlet_concentration_prior))
for k, dc in enumerate(dirichlet_concentration_prior):
    plot_results(plt.subplot(gs[0:2, k]), plt.subplot(gs[2, k]),
                 estimator, dc, X, y, plot_title=k == 0)

plt.show()
