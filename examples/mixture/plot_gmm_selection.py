"""
================================
Gaussian Mixture Model Selection
================================

This example shows that model selection can be performed with Gaussian Mixture
Models (GMM) using :ref:`information-theory criteria <aic_bic>`. Model selection
concerns both the covariance type and the number of components in the model.

In this case, both the Akaike Information Criterion (AIC) and the Bayes
Information Criterion (BIC) provide the right result, but we only demo the
latter as BIC is better suited to identify the true model among a set of
candidates. Unlike Bayesian procedures, such inferences are prior-free.

"""

# Authors: The scikit-learn developers
# SPDX-License-Identifier: BSD-3-Clause

# %%
# Data generation
# ---------------
#
# We generate two components (each one containing `n_samples`) by randomly
# sampling the standard normal distribution as returned by `numpy.random.randn`.
# One component is kept spherical yet shifted and re-scaled. The other one is
# deformed to have a more general covariance matrix.

import numpy as np

n_samples = 500
np.random.seed(0)
C = np.array([[0.0, -0.1], [1.7, 0.4]])
component_1 = np.dot(np.random.randn(n_samples, 2), C)  # general
component_2 = 0.7 * np.random.randn(n_samples, 2) + np.array([-4, 1])  # spherical

X = np.concatenate([component_1, component_2])

# %%
# We can visualize the different components:

import matplotlib.pyplot as plt

plt.scatter(component_1[:, 0], component_1[:, 1], s=0.8)
plt.scatter(component_2[:, 0], component_2[:, 1], s=0.8)
plt.title("Gaussian Mixture components")
plt.axis("equal")
plt.show()

# %%
# Model training and selection
# ----------------------------
#
# We vary the number of components from 1 to 6 and the type of covariance
# parameters to use:
#
# - `"full"`: each component has its own general covariance matrix.
# - `"tied"`: all components share the same general covariance matrix.
# - `"diag"`: each component has its own diagonal covariance matrix.
# - `"spherical"`: each component has its own single variance.
#

from sklearn.mixture import GaussianMixtureIC

gm_ic = GaussianMixtureIC(min_components=1, max_components=6, covariance_type="all")
gm_ic.fit(X)

# %%
# Plot the BIC scores
# -------------------
#
# To ease the plotting we can create a `pandas.DataFrame` from the results of
# the cross-validation done by the grid search. We re-inverse the sign of the
# BIC score to show the effect of minimizing it.

import pandas as pd

from sklearn.model_selection import ParameterGrid

param_grid = list(
    ParameterGrid(
        {
            "n_components": range(1, 7),
            "covariance_type": ["spherical", "tied", "diag", "full"],
        }
    )
)
df = pd.DataFrame(param_grid)
df.columns = ["Type of covariance", "Number of components"]
df["BIC score"] = gm_ic.criterion_
df.sort_values(by="BIC score").head()

# %%
import seaborn as sns

sns.catplot(
    data=df,
    kind="bar",
    x="Number of components",
    y="BIC score",
    hue="Type of covariance",
)
plt.show()

# %%
# In the present case, the model with 2 components and full covariance (which
# corresponds to the true generative model) has the lowest BIC score and is
# therefore selected by the grid search.
#
# Plot the best model
# -------------------
#
# We plot an ellipse to show each Gaussian component of the selected model. For
# such purpose, one needs to find the eigenvalues of the covariance matrices as
# returned by the `covariances_` attribute. The shape of such matrices depends
# on the `covariance_type`:
#
# - `"full"`: (`n_components`, `n_features`, `n_features`)
# - `"tied"`: (`n_features`, `n_features`)
# - `"diag"`: (`n_components`, `n_features`)
# - `"spherical"`: (`n_components`,)

from matplotlib.patches import Ellipse
from scipy import linalg

color_iter = sns.color_palette("tab10", 2)[::-1]
Y_ = gm_ic.predict(X)

fig, ax = plt.subplots()

for i, (mean, cov, color) in enumerate(
    zip(
        gm_ic.means_,
        gm_ic.covariances_,
        color_iter,
    )
):
    v, w = linalg.eigh(cov)
    if not np.any(Y_ == i):
        continue
    plt.scatter(X[Y_ == i, 0], X[Y_ == i, 1], 0.8, color=color)

    angle = np.arctan2(w[0][1], w[0][0])
    angle = 180.0 * angle / np.pi  # convert to degrees
    v = 2.0 * np.sqrt(2.0) * np.sqrt(v)
    ellipse = Ellipse(mean, v[0], v[1], angle=180.0 + angle, color=color)
    ellipse.set_clip_box(fig.bbox)
    ellipse.set_alpha(0.5)
    ax.add_artist(ellipse)

plt.title(
    f"Selected GMM: {gm_ic.covariance_type_} model, {gm_ic.n_components_} components"
)
plt.axis("equal")
plt.show()

from sklearn.metrics import adjusted_rand_score
from sklearn.mixture import GaussianMixture

# %%
# Comparison on a "double-cigar" dataset
# ---------------------------------------

# We now illustrate the behavior of
# :class:`~sklearn.mixture.GaussianMixtureIC` on a challenging
# anisotropic dataset consisting of two long, thin Gaussian
# components oriented at ±45° ("crossing double cigar"). In this
# configuration, EM with a single random initialization can
# converge to a poor partition, while the Mahalanobis–Ward
# hierarchical initialization used inside GaussianMixtureIC
# provides a more stable clustering. We quantify this with the
# Adjusted Rand Index (ARI) against the known ground truth.


def make_crossing_double_cigar(
    n_samples=600,
    sep=3.0,
    var_long=4.0,
    var_short=0.05,
    random_state=1,
):
    """Two long, thin Gaussians crossing at ±45 degrees.

    The first component is elongated along +45°, the second along
    -45°. The means are placed at (-sep/2, 0) and (sep/2, 0).
    """
    rng = np.random.RandomState(random_state)
    n1 = n_samples // 2
    n2 = n_samples - n1

    base_cov = np.array([[var_long, 0.0], [0.0, var_short]])

    def rotation(theta):
        c, s = np.cos(theta), np.sin(theta)
        return np.array([[c, -s], [s, c]])

    R1 = rotation(np.deg2rad(45.0))
    R2 = rotation(np.deg2rad(-45.0))

    cov1 = R1 @ base_cov @ R1.T
    cov2 = R2 @ base_cov @ R2.T

    mean1 = np.array([-sep / 2.0, 0.0])
    mean2 = np.array([sep / 2.0, 0.0])

    X1 = rng.multivariate_normal(mean1, cov1, size=n1)
    X2 = rng.multivariate_normal(mean2, cov2, size=n2)
    X = np.vstack([X1, X2])
    y = np.array([0] * n1 + [1] * n2)

    return X, y


def plot_selected_gmm(model, X, ax, title, ari):
    """Reuse the ellipse plotting style from the main example."""
    n_components = len(model.means_)
    color_iter = sns.color_palette("tab10", n_components)[::-1]

    Y_ = model.predict(X)
    for i, (mean, cov, color) in enumerate(
        zip(model.means_, model.covariances_, color_iter)
    ):
        if not np.any(Y_ == i):
            continue

        ax.scatter(X[Y_ == i, 0], X[Y_ == i, 1], 0.8, color=color)

        # same eigen-decomposition logic as in the original example
        v, w = linalg.eigh(cov)
        angle = np.arctan2(w[0][1], w[0][0])
        angle = 180.0 * angle / np.pi  # convert to degrees
        v = 2.0 * np.sqrt(2.0) * np.sqrt(v)

        ellipse = Ellipse(mean, v[0], v[1], angle=180.0 + angle, color=color)
        ellipse.set_clip_box(ax.figure.bbox)
        ellipse.set_alpha(0.5)
        ax.add_artist(ellipse)

    ax.set_title(f"{title}\n(ARI = {ari:.2f})")
    ax.set_xlabel("Feature 1")
    ax.set_ylabel("Feature 2")
    ax.axis("equal")


# Generate the crossing double-cigar data
X_dc, y_true = make_crossing_double_cigar(
    n_samples=600,
    sep=3.0,
    var_long=4.0,
    var_short=0.05,
    random_state=1,
)

# Plain GaussianMixture with a single random initialization
gm_plain = GaussianMixture(
    n_components=2,
    covariance_type="full",
    init_params="random",
    n_init=1,
    random_state=0,
)
gm_plain.fit(X_dc)
labels_plain = gm_plain.predict(X_dc)
ari_plain = adjusted_rand_score(y_true, labels_plain)

# GaussianMixtureIC uses Mahalanobis–Ward hierarchical initialization
# internally before running EM and selecting the best model by BIC.
gm_ic = GaussianMixtureIC(
    min_components=2,
    max_components=2,
    covariance_type="full",
    random_state=0,
)
labels_ic = gm_ic.fit_predict(X_dc)
ari_ic = adjusted_rand_score(y_true, labels_ic)

fig, axes = plt.subplots(1, 2, figsize=(10, 4))

plot_selected_gmm(
    gm_plain,
    X_dc,
    ax=axes[0],
    title="GaussianMixture",
    ari=ari_plain,
)

plot_selected_gmm(
    gm_ic,
    X_dc,
    ax=axes[1],
    title="GaussianMixtureIC",
    ari=ari_ic,
)

plt.tight_layout()
plt.show()
