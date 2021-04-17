"""
===============================================================
Generalized Additive Model Bagging with Gradient Boosting Trees
===============================================================

.. currentmodule:: sklearn

In this example, we show how a Generalized Additive Model
(`GAMBoostingRegressor`) can be used to learn from a dataset generated from a
summation of different shape functions.
"""
print(__doc__)

# %%
# Generate dataset as a summation of different shape functions
# ------------------------------------------------------------
import numpy as np

rng = np.random.RandomState(42)
n_samples = 10_000

X = rng.uniform(0.1, 4, size=(n_samples, 6))

functions = [
    lambda x: x,
    lambda x: x**2,
    lambda x: 3 * np.log(x),
    lambda x: 2 * np.sin(x)
]

y = np.zeros(n_samples)
for i, func in enumerate(functions):
    y += func(X[:, i])

y += rng.randn(n_samples)

# %%
# Train a GAMBoostingRegressor
# ----------------------------
from sklearn.ensemble import GAMBoostingRegressor
gam_reg = GAMBoostingRegressor(random_state=0)
gam_reg.fit(X, y)

# %%
# Compare the GAM predicts to actual shape functions
# --------------------------------------------------
import matplotlib.pyplot as plt

fig, axes = plt.subplots(2, 2, figsize=(12, 8))
for feature_idx, ax in enumerate(axes.ravel()):
    out = gam_reg.apply(X, feature_idx=feature_idx)
    X_sort_idx = X[:, feature_idx].argsort()
    X_sorted = X[:, feature_idx][X_sort_idx]
    ax.plot(X_sorted, out[X_sort_idx], label="predicted", alpha=0.5)
    ax.plot(X_sorted, functions[feature_idx](X_sorted), label="actual",
            alpha=0.5)
    ax.set_title(f"X{feature_idx}")
    ax.legend()

plt.show()
