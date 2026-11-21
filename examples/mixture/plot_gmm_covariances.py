"""
===============
GMM Covariances
===============

This example illustrates the impact of different covariance parameterizations
in :class:`~sklearn.mixture.GaussianMixture`.

The choice of covariance structure directly constrains the geometry of the
clusters and affects both model expressivity and generalization.
"""

# Authors: The scikit-learn developers
# SPDX-License-Identifier: BSD-3-Clause

# %%
# Covariance parameterizations
# ----------------------------
#
# Each covariance type defines how covariance matrices are shared or constrained
# across Gaussian components:
#
# - **full**: each component has its own general covariance matrix.
# - **tied**: all components share the same general covariance matrix.
# - **diag**: each component has its own diagonal covariance matrix.
# - **spherical**: each component has its own single variance.
#
# More flexible covariance models can capture more complex structure in the data,
# but they also require more parameters and can be less robust when data is limited.

# %%
# Dataset
# -------
#
# We construct a synthetic dataset drawn from a mixture of four Gaussians:
#
# - **Component 1**: spherical covariance
# - **Component 2**: diagonal covariance
# - **Component 3**: full covariance
# - **Component 4**: spherical covariance (fewer datapoints)

import numpy as np

rnd = np.random.RandomState(4)

n_components = 4
true_means = [(0, 0), (-2, 2), (-3, -2), (5, 0)]
true_covs = [np.eye(2), np.diag([1, 5]), [[3, 1], [1, 3]], np.eye(2)]
true_weights = [0.3, 0.3, 0.3, 0.1]

true_labels = rnd.choice(n_components, size=100, p=true_weights)

points = np.empty((true_labels.size, 2))

for label, mean, cov in zip(range(n_components), true_means, true_covs):
    idx = true_labels == label
    points[idx] = rnd.multivariate_normal(mean, cov, idx.sum())

# %%
# Model fitting
# --------------
#
# We fit four Gaussian mixture models with identical hyperparameters, varying only
# the covariance type.

from sklearn.mixture import GaussianMixture

gmms = {cov_type: {} for cov_type in ["spherical", "diag", "tied", "full"]}

for cov_type, gmm in gmms.items():
    model = GaussianMixture(
        n_components=n_components, covariance_type=cov_type, max_iter=20, random_state=0
    ).fit(points)

    gmm["model"] = model

# %%
# Label alignment and parameter extraction
# -----------------------------------------
#
# Since Gaussian mixture components are permutation-invariant, predicted labels
# must be aligned with ground truth for meaningful comparison.
#
# We also reconstruct the full covariance matrix from the packed representation.

from scipy.optimize import linear_sum_assignment

from sklearn.metrics.cluster import contingency_matrix

for cov_type, gmm in gmms.items():
    model = gmm["model"]

    labels = model.predict(points)
    cm = contingency_matrix(true_labels, labels)
    mapping = linear_sum_assignment(-cm)[1]

    gmm["labels"] = mapping[labels]
    gmm["means"] = model.means_[mapping]

    if cov_type == "full":
        covariances = model.covariances_

    elif cov_type == "tied":
        covariances = np.array([model.covariances_] * model.n_components)

    elif cov_type == "diag":
        covariances = np.array([np.diag(cov) for cov in model.covariances_])

    elif cov_type == "spherical":
        dim = model.means_.shape[1]
        covariances = np.array([np.eye(dim) * cov for cov in model.covariances_])

    gmm["covariances"] = covariances[mapping]

# %%
# Visualization of Gaussian components
# -------------------------------------
#
# Each Gaussian component is shown as an ellipse that summarizes the structure of
# the cluster.
#
# - **Center**: located at the mean of the component.
# - **Shape**: reflects how samples are distributed along different directions.
# - **Orientation**: indicates how the features vary together.
#
# If the ellipse is tilted, it means the variables are correlated.

import matplotlib as mpl
import matplotlib.pyplot as plt


def make_ellipse(mean, cov):
    v, w = np.linalg.eigh(cov)
    u = w[:, 0] / np.linalg.norm(w[:, 0])
    angle = np.degrees(np.arctan2(u[1], u[0]))
    v = 2.0 * np.sqrt(2.0) * np.sqrt(v)
    return mpl.patches.Ellipse(mean, v[0], v[1], angle=180 + angle)


def plot_gmm(ax, labels, means, covariances, colors):
    for label, (mean, cov) in enumerate(zip(means, covariances)):
        ell = make_ellipse(mean, cov)
        ell.set(alpha=0.5, color=colors[label], zorder=0)
        ax.add_patch(ell)
    ax.scatter(
        points[:, 0],
        points[:, 1],
        color=colors[labels],
        sizes=np.repeat(7, points.shape[0]),
        zorder=1,
    )


colors = np.array(["blue", "orange", "green", "magenta"])

fig, axes = plt.subplots(
    2, 3, subplot_kw={"xticks": [], "yticks": [], "aspect": "equal"}
)
axes = iter(axes.ravel())

ax0 = next(axes)
plot_gmm(ax0, true_labels, true_means, true_covs, colors)
ax0.set_title("Ground truth")

for (covariance_type, gmm), ax in zip(gmms.items(), axes):
    plot_gmm(ax, gmm["labels"], gmm["means"], gmm["covariances"], colors)
    ax.set_title(covariance_type)

for ax in axes:
    ax.set_visible(False)

fig.tight_layout()
plt.show()

# %%
# Evaluation
# ----------
#
# We evaluate each model using two complementary criteria:
#
# - **Rand index**: clustering agreement
# - **2-Wasserstein distance**: geometric distance between Gaussian components

import pandas as pd
from scipy.linalg import sqrtm

from sklearn.metrics import rand_score


def gaussian_wasserstein_distance(mean1, cov1, mean2, cov2):
    mean_term = np.sum((mean1 - mean2) ** 2)
    sqrt_cov1 = sqrtm(cov1)
    sqrt_middle = sqrtm(sqrt_cov1 @ cov2 @ sqrt_cov1)
    cov_term = np.trace(cov1 + cov2 - 2 * np.real(sqrt_middle))
    return np.sqrt(mean_term + max(cov_term, 0.0))


rows = {}

for covariance_type, gmm in gmms.items():
    row = {}

    row["RI"] = rand_score(true_labels, gmm["labels"])

    for true_mean, true_cov, pred_mean, pred_cov, color in zip(
        true_means, true_covs, gmm["means"], gmm["covariances"], colors
    ):
        w2 = gaussian_wasserstein_distance(true_mean, true_cov, pred_mean, pred_cov)
        row[f"W-{color}"] = w2

    rows[covariance_type] = row

df = pd.DataFrame.from_dict(rows, orient="index")

print(df.to_string(float_format="%.2f"))


# %%
# Interpretation
# ---------------
#
# Full covariance models are highly flexible and can
# accurately recover complex cluster shapes. However, they may overfit,
# especially for weakly represented clusters (e.g. magenta).
#
# In contrast, spherical covariance models impose isotropic constraints
# and ignore feature correlations,
# behaving similarly to
# `K-Means <https://scikit-learn.org/stable/auto_examples/cluster/plot_kmeans_assumptions.html>`_,
# and producing decision boundaries resembling to Voronoi partitions.
