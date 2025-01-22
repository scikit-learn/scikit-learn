# %%
import os

import numpy as np
import torch

import sklearn
from sklearn.datasets import make_blobs
from sklearn.mixture import GaussianMixture

os.environ["SCIPY_ARRAY_API"] = "1"

X, y = make_blobs(n_samples=int(1e3), n_features=2, centers=3, random_state=0)
X, y = torch.asarray(X), torch.asarray(y)

sklearn.set_config(array_api_dispatch=True)

gmm = GaussianMixture(
    n_components=3, covariance_type="full", random_state=0, init_params="random"
).fit(X)
print(gmm.means_)
print(gmm.covariances_)

# %%
import matplotlib as mpl
import matplotlib.pyplot as plt

fig, ax = plt.subplots()

ax.scatter(X[:, 0], X[:, 1], c=y)


def make_ellipses(gmm, ax):
    colors = ["navy", "turquoise", "darkorange"]
    for n, color in enumerate(colors):
        if gmm.covariance_type == "full":
            covariances = gmm.covariances_[n][:2, :2]
        elif gmm.covariance_type == "tied":
            covariances = gmm.covariances_[:2, :2]
        elif gmm.covariance_type == "diag":
            covariances = np.diag(gmm.covariances_[n][:2])
        elif gmm.covariance_type == "spherical":
            covariances = np.eye(gmm.means_.shape[1]) * gmm.covariances_[n]
        v, w = np.linalg.eigh(covariances)
        u = w[0] / np.linalg.norm(w[0])
        angle = np.arctan2(u[1], u[0])
        angle = 180 * angle / np.pi  # convert to degrees
        v = 2.0 * np.sqrt(2.0) * np.sqrt(v)
        ell = mpl.patches.Ellipse(
            gmm.means_[n, :2], v[0], v[1], angle=180 + angle, color=color
        )
        ell.set_clip_box(ax.bbox)
        ell.set_alpha(0.5)
        ax.add_artist(ell)
        ax.set_aspect("equal", "datalim")


make_ellipses(gmm, ax)
