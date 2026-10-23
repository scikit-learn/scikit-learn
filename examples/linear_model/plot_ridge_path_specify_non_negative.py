"""
=====================================================
Plot Ridge coefficients with non-negative constraints
=====================================================

Shows the effect of collinearity in the coefficients of an estimator with the influence
of non-negative constraints. This example is an extension of
:ref:`sphx_glr_auto_examples_linear_model_plot_ridge_path.py` that includes the
`positive` parameter.

Each curve is one coefficient of :class:`~sklearn.linear_model.Ridge` on the 10x10
Hilbert matrix, plotted against the regularization ``alpha``. The ``positive`` parameter
chooses which coefficients must be greater than or equal to zero.

``positive=False`` leaves every coefficient free. On this ill-conditioned matrix those
coefficients change sign and grow as ``alpha`` decreases.

A boolean array of shape ``(n_features,)`` applies the bound only where an entry is
``True``. A ``False`` entry leaves that coefficient free, so it can be negative.
``positive=True`` bounds every coefficient.

For the smallest values of ``alpha``, the constrained paths level off. Coefficients that
the bound keeps from becoming negative stay at zero. The unconstrained path keeps
changing.

"""

# Authors: The scikit-learn developers
# SPDX-License-Identifier: BSD-3-Clause

import matplotlib.pyplot as plt
import numpy as np

from sklearn.linear_model import Ridge

# X is the 10x10 Hilbert matrix
X = 1.0 / (np.arange(1, 11) + np.arange(0, 10)[:, np.newaxis])
y = np.ones(10)

# %%
# Non-negative constraints
# ------------------------
#
# ``True`` marks a coefficient that must be greater than or equal to zero.
# Here that is every other feature, starting with the first.

positive_mask = np.arange(X.shape[1]) % 2 == 0

# %%
# Compute paths
# -------------
#
# ``solver="auto"`` uses L-BFGS-B when at least one coefficient is
# constrained, and a closed-form solver otherwise.

n_alphas = 200
alphas = np.logspace(-10, 0, n_alphas)

constraints = [
    ("No non-negative constraint", False),
    ("Alternating non-negative constraint", positive_mask),
    ("All coefficients non-negative", True),
]

paths = []
for _, positive in constraints:
    coefs = [
        Ridge(alpha=alpha, fit_intercept=False, positive=positive).fit(X, y).coef_
        for alpha in alphas
    ]
    paths.append(np.asarray(coefs))

# %%
# Display results
# ---------------
#
# At ``alpha=1`` the three solutions nearly agree, because the penalty
# keeps every coefficient positive. As ``alpha`` decreases, free
# coefficients spread apart. Constrained coefficients that would cross
# below zero stay on the bound, and those two paths flatten out.

fig, axes = plt.subplots(1, 3, figsize=(12, 4), sharex=True, sharey=True)

for ax, (title, _), coefs in zip(axes, constraints, paths):
    ax.plot(alphas, coefs)
    ax.set_xscale("log")
    ax.set_xlim(ax.get_xlim()[::-1])  # reverse axis
    ax.set_xlabel("alpha")
    ax.set_title(title)
    ax.axhline(0.0, color="0.4", linewidth=0.8)
    ax.margins(x=0)

axes[0].set_ylabel("weights")
fig.legend(
    [f"Feature {i + 1}" for i in range(X.shape[1])],
    loc="lower center",
    ncol=5,
    bbox_to_anchor=(0.5, -0.08),
    fontsize="small",
)
fig.tight_layout()
plt.show()
