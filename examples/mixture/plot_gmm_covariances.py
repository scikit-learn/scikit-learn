"""
===============
GMM Covariances
===============

This example illustrates the impact of different covariance parameterizations
in :class:`~sklearn.mixture.GaussianMixture`.

The choice of covariance structure directly constrains the geometry of the
clusters and affects both model expressivity and generalization.
"""

# %%
# Authors: The scikit-learn developers
# SPDX-License-Identifier: BSD-3-Clause

# %%
# Overview
# --------
#
# We compare four covariance types:
#
# - **Spherical**: each component has covariance :math:`\sigma^2 I`
# - **Diagonal**: covariance matrices are diagonal
# - **Tied**: all components share the same covariance matrix
# - **Full**: each component has its own general covariance matrix
#
# While the *full* model is the most expressive, it introduces more parameters,
# which increases the risk of overfitting for small datasets

# %%
import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np
from scipy.optimize import linear_sum_assignment

from sklearn.metrics import confusion_matrix, rand_score
from sklearn.mixture import GaussianMixture

# %%
# Dataset
# ------
#
# We construct a synthetic dataset drawn from a mixture of four Gaussians:
#
# - Component 1: spherical covariance
# - Component 2: diagonal covariance
# - Component 3: full covariance
# - Component 4: spherical covariance (fewer datapoints)

rnd = np.random.RandomState(4)

ground_means = [(0, 0), (-2, 2), (-3, -2), (5, 0)]

ground_covs = [np.eye(2), np.diag([1, 5]), [[3, 1], [1, 3]], np.eye(2)]

ground_sample_sizes = [50, 50, 50, 10]

colors = ["blue", "orange", "green", "magenta"]

normals = [
    rnd.multivariate_normal(m, c, sz)
    for m, c, sz in zip(ground_means, ground_covs, ground_sample_sizes)
]

n_components = len(normals)

ground_labels = np.concatenate([[i] * sz for i, sz in enumerate(ground_sample_sizes)])

X = np.vstack(normals)

# %%
# Model Fitting
# -------------
#
# We fit four GMMs with identical hyperparameters but different covariance types.

cov_types = ["spherical", "diag", "tied", "full"]

estimators = [
    GaussianMixture(
        n_components=n_components, covariance_type=cov_type, max_iter=20, random_state=0
    )
    for cov_type in cov_types
]

gmm_labels = [est.fit_predict(X) for est in estimators]

# %%
# Visualization
# -------------

# %% Visualizing the distribution
# ----------------------------
#
# Each Gaussian component is represented as an ellipse that summarizes how the
# data is spread around its center. The center of the ellipse is the mean of the
# distribution. The shape of the ellipse comes from the covariance matrix: it
# stretches more in directions where the data varies a lot.


def make_ellipse(mean, cov):
    v, w = np.linalg.eigh(cov)
    u = w[:, 0] / np.linalg.norm(w[:, 0])
    angle = np.degrees(np.arctan2(u[1], u[0]))
    v = 2.0 * np.sqrt(2.0) * np.sqrt(v)
    return mpl.patches.Ellipse(mean, v[0], v[1], angle=180 + angle)


# %%
# Assigning correct labels
# ------------------------
#
# The GMM will cluster data into arbitrary clusters.
# We find the permutations to assign correct colors.


def make_perm(true_labels, pred_labels):
    cm = confusion_matrix(true_labels, pred_labels)
    return np.argsort(linear_sum_assignment(-cm)[1])


gmm_perms = [make_perm(ground_labels, gmm_cluster) for gmm_cluster in gmm_labels]

# %%
# Retrieving covariance matrices
# ------------------------------
#
# We convert all covariance formats to full matrices:


def get_covariances(gmm):
    if gmm.covariance_type == "full":
        return gmm.covariances_
    elif gmm.covariance_type == "tied":
        return [gmm.covariances_] * gmm.n_components
    elif gmm.covariance_type == "diag":
        return [np.diag(c) for c in gmm.covariances_]
    elif gmm.covariance_type == "spherical":
        return [np.eye(gmm.means_.shape[1]) * c for c in gmm.covariances_]


# %%
# Results
# -------
#
# We compare the learned cluster shapes against the ground truth.

fig = plt.figure(constrained_layout=True)
gs = fig.add_gridspec(2, 3)
gs.tight_layout(fig, pad=1.05)
axes = []
for i in range(len(estimators) + 1):
    axes.append(fig.add_subplot(gs[i // 3, i % 3]))

titles = ("Ground Truth", "Spherical", "Diagonal", "Tied", "Full")
all_means = [ground_means] + [gmm.means_ for gmm in estimators]
all_covs = [ground_covs] + [get_covariances(gmm) for gmm in estimators]
all_labels = [ground_labels] + gmm_labels
all_perms = [range(n_components)] + gmm_perms

for ax, title, means, covs, labels, perms in zip(
    axes, titles, all_means, all_covs, all_labels, all_perms
):
    ax.set(title=title, xticks=[], yticks=[], aspect="equal")
    score = rand_score(all_labels[0], labels)
    ax.text(
        0.05,
        0.95,
        f"RI = {score:.2f}",
        transform=ax.transAxes,
        va="top",
        ha="left",
        bbox=dict(boxstyle="round", facecolor="white", alpha=0.7),
    )

    clustered = [X[labels == k] for k in range(n_components)]

    for n, mean, cov, pts in zip(perms, means, covs, clustered):
        ell = make_ellipse(mean, cov)
        ell.set(alpha=0.5, color=colors[n], zorder=0)
        ax.add_patch(ell)
        ax.scatter(
            pts[:, 0],
            pts[:, 1],
            color=colors[n],
            sizes=np.repeat(7, pts.shape[1]),
            zorder=1,
        )

plt.show()

# %%
#
# As can be seen in the plots, each covariant type forces the shape of the
# cluster.
#
# Notice that the poor performance of the spherical type. It suffers similar
# pitfalls as :class:`~sklearn.cluster.KMeans`. However the predicted
# distribution for the sparse purple cluster matches the ground truth much more
# closely.
