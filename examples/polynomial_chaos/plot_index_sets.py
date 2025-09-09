r"""
===============================================================================
Polynomial Chaos regression
===============================================================================

This example illustrates how to use Polynomial Chaos regression for a simple
polynomial model given in [Saltelli2000]_. The model is defined as

.. math::
    y = \frac{1}{2^d}\prod_{j=1}^d (3x_j^2 + 1)

where the input variables :math:`x_j` follow a standard uniform distribution.

.. topic: References

    .. [Saltelli2000] `Saltelli, Andrea, et al. "Global sensitivity analysis:
       the primer." John Wiley & Sons, 2008.`

"""

# Authors: The scikit-learn developers
# SPDX-License-Identifier: BSD-3-Clause

# %%
# Let's fix the random seed for reproducibility.
import numpy as np

np.random.seed(2023)

# %%
# Consider a two-dimensional model (:math:`d = 2`).
dimension = 2

# %%
# Since the model is a product of :math:`d` two-dimensional polynomials, we
# should be able to exactly represent it using a polynomial of degree :math:`2
# \times d = 4`.
degree = 4

# %%
# The inputs :math:`x_1` and :math:`x_2` are assumed to follow a standard
# uniform distribution, i.e., :math:`x_j\sim\mathcal{U}[0, 1]` for :math:`j=1,
# 2`. First, let's generate some training data by sampling from those
# distributions.
from scipy.stats import uniform

distribution = uniform()
X = distribution.rvs((50, dimension))

# %%
# Then, evaluate the model in the training samples.
y = np.prod((3 * X**2 + 1) / 2, axis=1)

# %%
# Next, we fit a Polynomial Chaos surrogate to the available data.
from sklearn.polynomial_chaos import PolynomialChaosRegressor

pce = PolynomialChaosRegressor(distribution, degree=degree)
pce.fit(X, y)

# %%
# After the Polynomial Chaos expansion is fitted, we can visualize the terms
# in the basis, as well as the value of the coefficient associated with that
# basis term. The function below can be used to visualize the polynomial basis
# terms in two dimensions.
import matplotlib.pyplot as plt
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
    ax.set_xlabel("$x_1$")
    ax.yaxis.set_major_locator(MaxNLocator(integer=True))
    ax.set_ylabel("$x_2$")
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    ax.spines["bottom"].set_visible(False)
    ax.spines["left"].set_visible(False)
    ax.tick_params(length=0)


# %%
# Next, we plot the coefficients in our Polynomial Chaos expansion. Each square
# in the plot below represents a basis term. The coordinates of the square
# agrees with the degrees of the orthogonal polynomials that are associated
# with that basis term. In particular, a square at :math:`(2, 0)` agrees with
# the polynomial basis term :math:`p_2(x_1)p_0(x_2)`, where :math:`p_j(x)`
# represents the :math:`j`\ th Legendre polynomial. The color of the square
# indicates the magnitude of the coefficient in the expansion that is
# associated with this basis term.
_, ax = plt.subplots()
plot_coefficients(pce, ax)
ax.set_title("coefficients in PC expansion")
plt.show()

# %%
# The default choice is to add polynomial basis terms according to a
# `total_degree` scheme, meaning that only basis terms for which the *total*
# degree (i.e., the sum of the degrees of the individual polynomials) does not
# exceed a certain value. A multiindex set of this type is a simplex, as
# indicated by the pattern above. There are 4 different predefined polynomial
# multiindex set shapes implemented: `full_tensor`, `total_degree`,
# `hyperbolic_cross` and `Zaremba_cross`. Those multiindex set shapes are
# visualized below.
truncations = [
    "full_tensor",
    "total_degree",
    "hyperbolic_cross",
    "Zaremba_cross",
]
_, axes = plt.subplots(2, 2)
for ax, truncation in zip(axes.ravel(), truncations):
    pce = PolynomialChaosRegressor(distribution, degree=6, truncation=truncation)
    pce.fit(X, y)
    plot_coefficients(pce, ax, s=100)
    ax.set_title(truncation.capitalize().replace("_", " "))
plt.tight_layout()

# %%
# In practice, however, it is hard to know a priori which multiindex set shape
# (and which maximum polynomial degree) will give the best result for a given
# set of training points. We can use the
# :class:`~sklearn.model_selection.GridSearchCV` method to help us select the
# best combination.
from sklearn.model_selection import GridSearchCV

param_grid = [
    {
        "degree": [0, 1, 2, 3, 4, 5, 6],
        "truncation": [
            "full_tensor",
            "total_degree",
            "hyperbolic_cross",
            "Zaremba_cross",
        ],
    }
]
pceCV = GridSearchCV(
    PolynomialChaosRegressor(distribution),
    param_grid,
    scoring="neg_root_mean_squared_error",
)
pceCV.fit(X, y)
pceCV.best_params_

# %%
# The best multiindex set shape is a full tensor index set with degree 2. This
# is obvious when recalling the specific form of the model problem, and also
# from the the plot above: coefficients outside the rectangle from :math:`(0,
# 0)` to :math:`(2, 2)` are non-zero.

# %%
# See also
#   * :ref:`sphx_glr_auto_examples_polynomial_chaos_plot_pcr_ishigami.py` for
#     an example of how to extract information about the sensitivity of the
#     model output to the input features.
