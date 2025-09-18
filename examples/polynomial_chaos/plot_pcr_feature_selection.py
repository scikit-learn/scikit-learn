"""
===============================================================================
Pruning basis terms in Polynomial Chaos regression
===============================================================================

This example illustrates how to use Polynomial Chaos regression with basis
pruning during fitting. The pruning is performed using a feature selection
method, which discards irrelevant polynomial terms from the expansion. This
can lead to more interpretable models and reduces computational cost during
prediction.

"""

# Authors: The scikit-learn developers
# SPDX-License-Identifier: BSD-3-Clause

# %%
# Let's fix the random seed for reproducibility.
import numpy as np

np.random.seed(1234)

# %%
# We consider the following polynomial model with interaction:
#
# .. math::
#   f(x_1, x_2) = 3 x_1^2 + 2 x_2 + 0.5 x_1 x_2.
#
# We also add a small amount of Gaussian noise to the outputs.
X = np.random.rand(100, 2)
y = 3 * X[:, 0] ** 2 + 2 * X[:, 1] + 0.5 * X[:, 0] * X[:, 1]
y += 0.1 * np.random.randn(100)

# %%
# The inputs :math:`x_1, x_2` are assumed to follow a uniform distribution
# on [0, 1].
from scipy.stats import uniform

distribution = [uniform()] * 2

# %%
# Let's split the data into training and test sets.
from sklearn.model_selection import train_test_split

X_train, X_test, y_train, y_test = train_test_split(X, y, random_state=0)

# %%
# We first fit a Polynomial Chaos expansion of degree 5 without pruning.
from sklearn.linear_model import LinearRegression
from sklearn.polynomial_chaos import PolynomialChaosRegressor

pce = PolynomialChaosRegressor(
    distribution=distribution, degree=5, solver=LinearRegression(fit_intercept=False)
)
pce.fit(X_train, y_train)
print(f"R^2 on test set: {pce.score(X_test, y_test):.2f}")

# %%
# Let's plot the learned coefficients.
import matplotlib.pyplot as plt

indices = [str(m) for m in pce.multiindices_]
plt.figure(figsize=(8, 4))
plt.stem(range(len(pce.coef_)), pce.coef_)
plt.xticks(range(len(pce.coef_)), indices, rotation=90)
plt.xlabel("Polynomial basis multi-index")
plt.ylabel("Coefficient value")
plt.title("Polynomial Chaos coefficients (no pruning)")
plt.tight_layout()
plt.show()

# %%
# We see that many higher-order terms have nonzero coefficients, even though
# the true model only depends on a few terms.
#
# To remove irrelevant basis terms, we can use a feature selector. Here, we
# use :class:`~sklearn.feature_selection.SelectFromModel` with
# :class:`~sklearn.linear_model.LassoCV` as the underlying estimator.
from sklearn.feature_selection import SelectFromModel
from sklearn.linear_model import LassoCV

pruned_pce = PolynomialChaosRegressor(
    distribution=distribution,
    degree=5,
    solver=LinearRegression(fit_intercept=False),
    feature_selector=SelectFromModel(LassoCV(fit_intercept=False), max_features=5),
)
pruned_pce.fit(X_train, y_train)
print(f"R^2 on test set: {pruned_pce.score(X_test, y_test):.2f}")

# %%
# Let's inspect the selected terms. The rows correspond to multi-indices
# :math:`(\alpha_1, \alpha_2)`, indicating which powers of :math:`x_1` and
# :math:`x_2` are included in the expansion.
print("Selected terms:")
for mindex in pruned_pce.multiindices_:
    print(mindex)

# %%
# Pruning recovers the correct sparse structure of the model
# (quadratic in :math:`x_1`, linear in :math:`x_2`, and with an interaction
# term), while maintaining a high predictive accuracy on the test data.
from matplotlib.cm import ScalarMappable
from matplotlib.colors import Normalize
from matplotlib.ticker import MaxNLocator


def plot_coefficients(pce, ax, s=600):
    cmap = plt.get_cmap("Spectral_r")
    cmi = min(np.abs(pce.coef_.flatten()))
    cma = max(np.abs(pce.coef_.flatten()))
    for idx, c in zip(pce.multiindices_, np.abs(pce.coef_.flatten())):
        ax.scatter(*idx, color=cmap((c - cmi) / (cma - cmi)), marker="s", s=s)
    norm = Normalize(vmin=cmi, vmax=cma)
    plt.colorbar(
        ScalarMappable(norm=norm, cmap=cmap),
        ax=ax,
        label="value of the coefficient",
    )
    ax.xaxis.set_major_locator(MaxNLocator(integer=True))
    ax.set_xlabel("power of $x_1$")
    ax.yaxis.set_major_locator(MaxNLocator(integer=True))
    ax.set_ylabel("power of $x_2$")
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    ax.spines["bottom"].set_visible(False)
    ax.spines["left"].set_visible(False)
    ax.tick_params(length=0)


# %%
# Let's plot the learned coefficients in multi-index space.
fig, ax = plt.subplots(figsize=(5, 4))
plot_coefficients(pce, ax)
ax.set_title("Polynomial Chaos coefficients (no pruning)")
plt.show()

# %%
# After pruning:
fig, ax = plt.subplots(figsize=(5, 4))
plot_coefficients(pruned_pce, ax)
ax.set_title("Polynomial Chaos coefficients (with pruning)")
ax.set_xlim(-0.5, 5.5)
ax.set_ylim(-0.5, 5.5)
plt.show()

# %%
# See also
#   * :ref:`sphx_glr_auto_examples_polynomial_chaos_plot_pcr_noisy.py` for
#     an example on handling noisy outputs.
#   * :ref:`sphx_glr_auto_examples_polynomial_chaos_plot_pcr_sobol_g.py` for
#     an example of adaptive construction of multiindex sets.
