# %%

"""
============================================================
Custom initialization for GaussianMixture with Ward-Mahalanobis
============================================================

This example illustrates how the EM algorithm for
:class:`~sklearn.mixture.GaussianMixture` can converge to different local optima
depending on initialization, and how to provide custom initial parameters.

We compare a built-in initialization strategy (k-means++) to a Ward-Mahalanobis
hierarchical initialization:

- Estimate a pooled covariance matrix on the data.
- Whiten the data so that Euclidean distances in the transformed space
  correspond to Mahalanobis distances in the original space.
- Run Ward clustering on the whitened data and cut it to obtain cluster labels.
- Turn those clusters into ``weights_init``, ``means_init``, and
  ``precisions_init`` for :class:`~sklearn.mixture.GaussianMixture`.
"""

# Authors: The scikit-learn developers
# SPDX-License-Identifier: BSD-3-Clause

# %%
# Imports
# -------

import matplotlib.pyplot as plt
import numpy as np
from matplotlib.patches import Ellipse
from scipy import linalg

from sklearn.cluster import AgglomerativeClustering
from sklearn.metrics import adjusted_rand_score
from sklearn.mixture import GaussianMixture

# %%
# Generate a dataset where Mahalanobis-aware initialization helps
# --------------------------------------------------------------
#
# We build two spherical Gaussians in a latent space, separated along one axis,
# and apply a strong linear transform so that separation happens mostly along a
# low-variance direction in the observed space.


def make_low_variance_separation(
    n_samples=400,
    sep=4.0,
    random_state=0,
):
    """Two spherical clusters separated along y, then strongly transformed."""
    rng = np.random.RandomState(random_state)
    n1 = n_samples // 2
    n2 = n_samples - n1

    y1 = rng.randn(n1, 2) + np.array([0.0, -sep / 2.0])
    y2 = rng.randn(n2, 2) + np.array([0.0, sep / 2.0])
    Y = np.vstack([y1, y2])
    y_true = np.array([0] * n1 + [1] * n2)

    # Strong anisotropy + correlation in observed space.
    A = np.array([[8.0, 2.0], [0.0, 0.08]])
    X = Y @ A.T

    return X, y_true


X, y_true = make_low_variance_separation(n_samples=400, sep=4.0, random_state=0)


# %%
# Whitening and Ward-Mahalanobis initialization
# --------------------------------------------


def whiten_with_pooled_covariance(X, reg_covar=1e-6):
    """Whiten centered X using the pooled covariance (Cholesky-based)."""
    X = np.asarray(X)
    X_mean = X.mean(axis=0)
    Xc = X - X_mean
    n_features = X.shape[1]

    cov = np.cov(Xc, rowvar=False)
    cov = np.atleast_2d(cov)
    cov.flat[:: n_features + 1] += reg_covar

    L = linalg.cholesky(cov, lower=True)  # cov = L L^T
    X_white = linalg.solve_triangular(L, Xc.T, lower=True).T

    return X_white


def mahalanobis_ward_initialization(
    X,
    n_components,
    reg_covar=1e-6,
):
    """Ward clustering in whitened space -> (weights, means, precisions)."""
    X = np.asarray(X)
    n_samples, n_features = X.shape

    X_white = whiten_with_pooled_covariance(X, reg_covar=reg_covar)
    labels = AgglomerativeClustering(
        n_clusters=n_components,
        linkage="ward",
    ).fit_predict(X_white)

    weights = np.bincount(labels, minlength=n_components).astype(float)
    weights /= float(n_samples)

    means = np.zeros((n_components, n_features), dtype=float)
    covariances = np.zeros((n_components, n_features, n_features), dtype=float)

    x_mean = X.mean(axis=0)
    Xc = X - x_mean
    global_cov = np.cov(Xc, rowvar=False)
    global_cov = np.atleast_2d(global_cov)
    global_cov.flat[:: n_features + 1] += reg_covar

    for k in range(n_components):
        Xk = X[labels == k]
        if Xk.shape[0] <= 1:
            means[k] = x_mean if Xk.shape[0] == 0 else Xk[0]
            Ck = global_cov
        else:
            means[k] = Xk.mean(axis=0)
            Ck = np.cov(Xk, rowvar=False)
            Ck = np.atleast_2d(Ck)
            Ck.flat[:: n_features + 1] += reg_covar
        covariances[k] = Ck

    precisions_init = np.empty_like(covariances)
    for k in range(n_components):
        precisions_init[k] = linalg.pinvh(covariances[k])

    return weights, means, precisions_init


# %%
# Visualize geometry: observed vs whitened
# ---------------------------------------

X_white = whiten_with_pooled_covariance(X, reg_covar=1e-6)

fig, axes = plt.subplots(1, 2, figsize=(10, 4), sharey=False)

for label, color in zip([0, 1], plt.cm.tab10([0, 1])):
    mask = y_true == label
    axes[0].scatter(X[mask, 0], X[mask, 1], s=3, color=color, label=f"Class {label}")
axes[0].set_title("Observed space (colored by true labels)")
axes[0].set_xlabel("Feature 1")
axes[0].set_ylabel("Feature 2")
axes[0].legend(loc="best", markerscale=3, frameon=False)

for label, color in zip([0, 1], plt.cm.tab10([0, 1])):
    mask = y_true == label
    axes[1].scatter(X_white[mask, 0], X_white[mask, 1], s=3, color=color)
axes[1].set_title("Whitened space (Ward-Mahalanobis)")
axes[1].set_xlabel("Whitened feature 1")
axes[1].set_ylabel("Whitened feature 2")
axes[1].set_aspect("equal", adjustable="box")

fig.tight_layout()
plt.show()


# %%
# Precompute custom init once
# ---------------------------

weights_init, means_init, precisions_init = mahalanobis_ward_initialization(
    X,
    n_components=2,
    reg_covar=1e-6,
)


# %%
# Compare initialization strategies (ARI across seeds)
# ---------------------------------------------------

seeds = range(10)
methods = [
    ("random (1 start)", {"init_params": "random", "n_init": 1}),
    ("random (5 starts)", {"init_params": "random", "n_init": 5}),
    ("kmeans", {"init_params": "kmeans", "n_init": 1}),
    ("k-means++", {"init_params": "k-means++", "n_init": 1}),
    ("random_from_data", {"init_params": "random_from_data", "n_init": 1}),
    ("Ward-Mahalanobis (custom)", None),
]
method_names = [name for name, _ in methods]

rows = []
for seed in seeds:
    for name, params in methods:
        if params is None:
            # Passing explicit init parameters makes EM deterministic here.
            gmm = GaussianMixture(
                n_components=2,
                covariance_type="full",
                reg_covar=1e-6,
                weights_init=weights_init,
                means_init=means_init,
                precisions_init=precisions_init,
            )
        else:
            gmm = GaussianMixture(
                n_components=2,
                covariance_type="full",
                reg_covar=1e-6,
                random_state=seed,
                **params,
            )

        gmm.fit(X)
        labels = gmm.predict(X)

        rows.append(
            {
                "Initialization": name,
                "Seed": seed,
                "ARI": adjusted_rand_score(y_true, labels),
            }
        )

# Convert to a dict-of-lists for simple plotting without pandas.
results_by_method = {name: [] for name in method_names}
for row in rows:
    results_by_method[row["Initialization"]].append(row["ARI"])

# Matplotlib jitter plot to show all points (and the median).
rng = np.random.RandomState(0)
fig, ax = plt.subplots(figsize=(10, 3.4))

for i, name in enumerate(method_names):
    y = np.asarray(results_by_method[name])
    x = i + 0.15 * (rng.rand(y.size) - 0.5)  # small horizontal jitter
    ax.scatter(x, y, s=25, alpha=0.8)

    med = float(np.median(y))
    ax.plot([i - 0.22, i + 0.22], [med, med], linewidth=3)

ax.set_title("Clustering quality across seeds (ARI)")
ax.set_ylabel("Adjusted Rand Index")
ax.set_ylim(-0.05, 1.05)
ax.set_xticks(range(len(method_names)))
ax.set_xticklabels(method_names, rotation=20, ha="right")
fig.tight_layout()
plt.show()


# %%
# Visualize a difficult run
# -------------------------
#
# We pick the seed giving the worst ARI for k-means++ and compare against the
# Ward-Mahalanobis custom initialization.
#
# The plotting style below mirrors the ellipse visualization used in
# ``plot_gmm_selection.py`` (adapted here for the fixed ``covariance_type="full"``
# setting).


def plot_gmm(ax, gmm, X, title):
    """Plot points colored by the GMM prediction and draw covariance ellipses."""
    Y = gmm.predict(X)
    colors = plt.cm.tab10(np.linspace(0, 1, gmm.n_components))[::-1]

    for i, (mean, color) in enumerate(zip(gmm.means_, colors)):
        if not np.any(Y == i):
            continue

        ax.scatter(X[Y == i, 0], X[Y == i, 1], s=3, color=color)

        cov = gmm.covariances_[i]
        v, w = linalg.eigh(cov)

        angle = np.arctan2(w[0, 1], w[0, 0])
        angle = 180.0 * angle / np.pi
        v = 2.0 * np.sqrt(2.0) * np.sqrt(v)

        ellipse = Ellipse(mean, v[0], v[1], angle=180.0 + angle, color=color, alpha=0.5)
        ellipse.set_clip_box(ax.figure.bbox)
        ax.add_artist(ellipse)

    ax.set_title(title)
    ax.set_xlabel("Feature 1")
    ax.set_ylabel("Feature 2")
    ax.set_aspect("auto")


# Find worst seed for k-means++ (using the rows list computed above).
kpp_aris = [(r["Seed"], r["ARI"]) for r in rows if r["Initialization"] == "k-means++"]
worst_seed = min(kpp_aris, key=lambda t: t[1])[0]

gmm_kpp = GaussianMixture(
    n_components=2,
    covariance_type="full",
    init_params="k-means++",
    n_init=1,
    reg_covar=1e-6,
    random_state=worst_seed,
).fit(X)
ari_kpp = adjusted_rand_score(y_true, gmm_kpp.predict(X))

gmm_ward = GaussianMixture(
    n_components=2,
    covariance_type="full",
    reg_covar=1e-6,
    weights_init=weights_init,
    means_init=means_init,
    precisions_init=precisions_init,
).fit(X)
ari_ward = adjusted_rand_score(y_true, gmm_ward.predict(X))

fig, axes = plt.subplots(1, 2, figsize=(10, 4), sharex=True, sharey=True)

plot_gmm(
    axes[0],
    gmm_kpp,
    X,
    title=f"GaussianMixture\n(k-means++, 1 start)\n(ARI={ari_kpp:.2f})",
)
plot_gmm(
    axes[1],
    gmm_ward,
    X,
    title=f"GaussianMixture\n(Ward-Mahalanobis init)\n(ARI = {ari_ward:.2f})",
)

fig.tight_layout()
plt.show()
