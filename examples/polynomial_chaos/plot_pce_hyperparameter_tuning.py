r"""
===============================================================================
Hyperparameter tuning of a Polynomial Chaos Expansion with GridSearchCV
===============================================================================

This example shows how to use :class:`~sklearn.model_selection.GridSearchCV` to
tune the hyperparameters of a Polynomial Chaos Expansion (PCE). We demonstrate
how cross-validation can identify an appropriate polynomial degree when the
data contains noise.
"""

# Authors: The scikit-learn developers
# SPDX-License-Identifier: BSD-3-Clause

# %%
# Reproducibility
import numpy as np
from scipy.stats import uniform

from sklearn.utils import check_random_state

random_state = check_random_state(123)

# %%
# Generate noisy data
# -------------------
# We define a nonlinear model:
#
# .. math::
#   y = 2 + 0.5 \sin(\pi x_1) + 0.3 x_2^2 + \epsilon,
#
# where :math:`x_1, x_2 \sim \mathcal{U}(-1, 1)` and
# :math:`\epsilon \sim \mathcal{N}(0, 0.05^2)`.
#
# The signal contains both polynomial and non-polynomial structure. With noise
# present, very high polynomial degrees will overfit.
X = uniform(loc=-1, scale=2).rvs((200, 2), random_state=random_state)
noise = random_state.normal(0, 0.05, size=len(X))
y = 2 + 0.5 * np.sin(np.pi * X[:, 0]) + 0.3 * X[:, 1] ** 2 + noise

# %%
# Define grid of hyperparameters
from sklearn.model_selection import GridSearchCV
from sklearn.polynomial_chaos import PolynomialChaosExpansion

param_grid = [
    {
        "degree": list(range(0, 9)),  # try degrees 0 to 8
        "truncation": ["total_degree"],
    }
]

# %%
# Run grid search
pceCV = GridSearchCV(
    PolynomialChaosExpansion(uniform(loc=-1, scale=2)),
    param_grid,
    scoring="neg_root_mean_squared_error",
    cv=5,
)
pceCV.fit(X, y)

print("Best parameters found:")
print(pceCV.best_params_)

# %%
# Inspect results
import pandas as pd

results = pd.DataFrame(pceCV.cv_results_)
cols = ["param_degree", "mean_test_score"]
results[cols].sort_values("mean_test_score", ascending=False).head()

# %%
# Visualization
# -------------
# Let's visualize the cross-validation RMSE as a function of degree,
# separately for LinearRegression and LassoCV.
import matplotlib.pyplot as plt

fig, ax = plt.subplots()
mean_scores = results.groupby("param_degree")["mean_test_score"].max()
ax.plot(mean_scores.index, -mean_scores.values, marker="o")

ax.set_xlabel("Polynomial degree")
ax.set_ylabel("Cross-validated RMSE")
ax.set_title("Grid search results for PCE")
ax.legend(frameon=False)
plt.show()

# %%
# See also
#   * :ref:`sphx_glr_auto_examples_polynomial_chaos_plot_pce_noisy_data.py` for
#     an example of how to use sparse estimators to deal with noisy measurements.
#   * :ref:`sphx_glr_auto_examples_polynomial_chaos_plot_pce_feature_selection_g.py`
#     for an example of how to use pruning to remove small basis terms from
#     the expansion.
# %%
