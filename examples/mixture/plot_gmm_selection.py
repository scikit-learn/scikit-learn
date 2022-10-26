"""
================================
Gaussian Mixture Model Selection
================================

This example shows that model selection can be performed with Gaussian Mixture
Models using :ref:`information-theory criteria <aic_bic>`. Model selection
concerns both the covariance type and the number of components in the model.

In this case, both the Akaike Information Criterion (AIC) and the Bayes
Information Criterion (BIC) provide the right result, but we only demo the
latter as BIC is better suited to identify the true model among a set of
candidates. Unlike Bayesian procedures, such inferences are prior-free.

"""

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

X = np.r_[
    np.dot(np.random.randn(n_samples, 2), C),
    0.7 * np.random.randn(n_samples, 2) + np.array([-6, 3]),
]

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
# We score the different models and keep the best model (the lowest BIC).

from sklearn.mixture import GaussianMixture

lowest_bic = np.infty
bic = []
n_components_range = range(1, 7)
cov_types = ["spherical", "tied", "diag", "full"]
for cov_type in cov_types:
    for n_components in n_components_range:
        # Fit a Gaussian mixture with EM
        gmm = GaussianMixture(n_components=n_components, covariance_type=cov_type)
        gmm.fit(X)
        bic.append(gmm.bic(X))
        if bic[-1] < lowest_bic:
            lowest_bic = bic[-1]
            best_gmm = gmm

bic = np.array(bic)
clf = best_gmm

# %%
# Plot the BIC scores
# -------------------
#
# We draw a `*` to highlight the model with the lowest BIC score. In the present
# case, the model with 2 components and full covariance (which corresponds to
# the true generative model) is selected.

import itertools
import matplotlib.pyplot as plt

color_iter = itertools.cycle(["navy", "darkorange", "cornflowerblue", "turquoise"])
bars = []

for i, (cov_type, color) in enumerate(zip(cov_types, color_iter)):
    xpos = np.array(n_components_range) + 0.2 * (i - 2)
    bars.append(
        plt.bar(
            xpos,
            bic[i * len(n_components_range) : (i + 1) * len(n_components_range)],
            width=0.2,
            color=color,
        )
    )
plt.xticks(n_components_range)
plt.ylim([bic.min() * 1.01 - 0.01 * bic.max(), bic.max()])
plt.title("BIC score per model")
xpos = (
    np.mod(bic.argmin(), len(n_components_range))
    + 0.65
    + 0.2 * np.floor(bic.argmin() / len(n_components_range))
)
plt.text(xpos, bic.min() * 0.97 + 0.03 * bic.max(), "*", fontsize=14)
plt.xlabel("Number of components")
plt.legend([b[0] for b in bars], cov_types)
plt.show()

# %%
# Plot the best model
# -------------------
#
# We plot an ellipse to show each Gaussian component of the selected model. For
# such purpose, one needs to find the eigenvalues of the covariance matrices as
# returned by the `covariances_` attribute. The shape of such matrices depends
# on the `covariance_type`:
#
# - `"full"`: (n_components, n_features, n_features)
# - `"tied"`: (n_features, n_features)
# - `"diag"`: (n_components, n_features)
# - `"spherical"`: (n_components,)

from matplotlib.patches import Ellipse
from scipy import linalg

fig, ax = plt.subplots()
Y_ = clf.predict(X)
for i, (mean, cov, color) in enumerate(zip(clf.means_, clf.covariances_, color_iter)):
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
    f"Selected GMM: {best_gmm.covariance_type} model, "
    f"{best_gmm.n_components} components"
)
plt.show()
