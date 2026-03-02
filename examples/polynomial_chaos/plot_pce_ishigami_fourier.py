r"""
===============================================================================
Exact PCE of the Ishigami function with a mixed Fourier-Legendre basis
===============================================================================

We revisit the Ishigami function from
:ref:`sphx_glr_auto_examples_polynomial_chaos_plot_pce_ishigami.py`

.. math::
    y = \sin x_1 + a \sin^2 x_2 + b \, x_3^4 \sin x_1,

with independent inputs uniformly distributed on :math:`[-\pi,\pi]`, and
parameters :math:`a=7`, :math:`b=0.1`.

This time we build a Polynomial Chaos Expansion with a *Fourier* basis in
:math:`x_1` and :math:`x_2` and a *Legendre* basis in :math:`x_3`. The
orthonormal Fourier system is

.. math::
    \psi_0(x)=1,\qquad
    \psi_{2k-1}(x)=\sqrt{2}\sin(kx),\qquad
    \psi_{2k}(x)=\sqrt{2}\cos(kx),

The model has **exactly five** nonzero product-basis terms. This example shows
how we can recover the *exact* sparse structure of the Ishigami function
automatically using a sparse estimator (such as
:class:`~sklearn.linear_model.LassoLarsCV`).
"""

# Authors: The scikit-learn developers
# SPDX-License-Identifier: BSD-3-Clause

# %%
# First, let's fix the random seed for reproducibility.
import numpy as np

np.random.seed(2025)

# %%
# Next, we define the Ishigami function.
#
# This is the same function as in the standard example, but here we will
# represent it in a mixed Fourier–Legendre basis. Recall that `a` and `b`
# are parameters that control the nonlinearity and interaction terms.
a = 7.0
b = 0.1


def ishigami(X):
    return (
        np.sin(X[:, 0]) + a * np.sin(X[:, 1]) ** 2 + b * X[:, 2] ** 4 * np.sin(X[:, 0])
    )


# %%
# Now, let's define the custom Fourier basis for variables `x1` and `x2`.
# This basis is orthonormal under the uniform distribution on [-pi, pi].
#
# We also wrap the distribution to give it a recognizable name ("fourier"),
# so it can be passed to the PolynomialChaosExpansion class alongside
# Legendre polynomials.
from scipy.stats import uniform

from sklearn.polynomial_chaos import PolynomialChaosExpansion
from sklearn.utils._orthogonal_polynomial import Polynomial


class Fourier(Polynomial):
    """Fourier basis orthonormal on [-pi, pi] under the uniform pdf."""

    def _vandermonde(self, points, degree):
        points = np.asarray(points)
        n = len(points)
        V = np.zeros((n, degree + 1))
        V[:, 0] = 1.0
        for j in range(1, degree + 1):
            k = (j + 1) // 2
            if j % 2 == 1:  # 1,3,5,... → sin(kx)
                V[:, j] = np.sqrt(2) * np.sin(k * points)
            else:  # 2,4,6,... → cos(kx)
                V[:, j] = np.sqrt(2) * np.cos(k * points)
        return V

    @staticmethod
    def _distribution():
        return "fourier"

    def scale_features_from_distribution(self, X, distribution):
        # Map Uniform(0,1) → [-pi, pi] so orthogonality is respected
        X = super().scale_features_from_distribution(X, distribution)
        return 2 * np.pi * (X - 0.5)

    def _norm_squared(self, degree):
        return 1  # orthonormal with respect to Uniform(-pi, pi)


class FourierDistribution:
    """A thin wrapper giving scipy Uniform a 'fourier' name."""

    def __init__(self, loc=-np.pi, scale=2 * np.pi):
        self._frozen = uniform(loc=loc, scale=scale)
        self.args = self._frozen.args
        self.kwds = self._frozen.kwds
        self.dist = type(
            "dist",
            (),
            {"name": "fourier", "_parse_args": self._frozen.dist._parse_args},
        )()

    def __getattr__(self, name):
        return getattr(self._frozen, name)


# %%
# With this, we can now define the joint distribution of our inputs:
#   - Fourier-distributed x1
#   - Fourier-distributed x2
#   - Uniform-distributed x3 (Legendre basis is orthonormal here)
U = uniform(loc=-np.pi, scale=2 * np.pi)
distribution = (
    FourierDistribution(),
    FourierDistribution(),
    U,
)

# %%
# Let's generate training data from the input distribution and evaluate
# the Ishigami function. We'll use 1000 samples for training.
X = U.rvs((1000, 3))
y = ishigami(X)

# %%
# Next, we fit a Polynomial Chaos Expansion (PCE) using the mixed basis.
# We allow terms up to total degree 5, but rely on the sparse regression
# estimator (LassoLarsCV) to automatically select only the significant ones.
#
# This approach will recover exactly the five active terms of the Ishigami
# function in this basis.
from sklearn.linear_model import LassoLarsCV

pce = PolynomialChaosExpansion(
    distribution=distribution,
    degree=5,
    scale_outputs=False,
    estimator=LassoLarsCV(fit_intercept=False),
)
pce.fit(X, y)

# %%
# To check accuracy, we compute the R^2 score on fresh test data.
# A value very close to 1 indicates that the PCE surrogate is essentially exact.
X_test = U.rvs((1024, 3))
y_test = ishigami(X_test)
print("R^2:", pce.score(X_test, y_test))

# %%
# Finally, let's inspect the recovered coefficients.
#
# Each basis function in the expansion is identified by a *multiindex*,
# for example:
#
#   [i, j, k]
#
# means:
#   - use the i-th univariate polynomial for `x1` (Fourier basis),
#   - the j-th for `x2` (Fourier basis),
#   - the k-th for `x3` (Legendre basis).
#
# The Fourier basis is
#
# .. math::
#   \psi_0(x) &= 1, \\
#   \psi_1(x) &= \sqrt{2}\,\sin(x), \\
#   \psi_2(x) &= \sqrt{2}\,\cos(x), \\
#   \psi_3(x) &= \sqrt{2}\,\sin(2x), \\
#   \psi_4(x) &= \sqrt{2}\,\cos(2x), \quad \ldots
#
# Therefore, the expected active multiindices are:
#
#   [0, 0, 0] → constant
#   [0, 4, 0] → cos(2·x2)
#   [1, 0, 0] → sin(x1)
#   [1, 0, 2] → sin(x1) * P2(x3)
#   [1, 0, 4] → sin(x1) * P4(x3)
#
# All other coefficients should vanish up to numerical tolerance.
for m, c in zip(pce.multiindices_, pce.coef_):
    if np.abs(c) > 1e-14:
        print(f"{m}: {c}")

# %%
# See also
#   * :ref:`sphx_glr_auto_examples_polynomial_chaos_plot_pce_ishigami.py` for
#     a global sensitivity analysis of the Ishigami function.
