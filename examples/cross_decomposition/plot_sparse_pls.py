"""
===================================
Sparse Partial Least Squares (SPLS)
===================================

This example demonstrates how to use Sparse Partial Least Squares (SPLS) with
a grid search for penalty parameters. SPLS is a dimensionality reduction
technique that combines the features of both Partial Least Squares (PLS) and
Lasso regression, allowing for sparsity in the resulting weight vectors.

We will compare SPLS with PLSCanonical to illustrate how SPLS can find true
sparsity in the weights.

"""

import matplotlib.pyplot as plt
import numpy as np

from sklearn.cross_decomposition import SPLS, PLSCanonical
from sklearn.model_selection import GridSearchCV

# Setting up a random state for reproducibility
n = 500
rng = np.random.RandomState(42)

# Generating two latent variables to represent the true underlying signal
l1 = rng.normal(size=n)
l2 = rng.normal(size=n)

# Combining the latent variables to form a matrix
latents = np.array([l1, l2]).T

# Generating additional noise columns to simulate non-informative features
noise_columns_X = rng.normal(size=4 * n).reshape((n, 4))
noise_columns_Y = rng.normal(size=4 * n).reshape((n, 4))

# Concatenating the latent variables with noise columns to form datasets X and Y
# This helps in testing the SPLS model's ability to identify the true signal
X = np.concatenate([latents, noise_columns_X], axis=1)
Y = np.concatenate([latents, noise_columns_Y], axis=1)

# Splitting the data into training and test sets
X_train = X[: n // 2]
Y_train = Y[: n // 2]
X_test = X[n // 2 :]
Y_test = Y[n // 2 :]

# Creating SPLS and PLSCanonical models
# SPLS model uses a grid search to find the best penalty parameters,
# illustrating its ability to induce sparsity in the model.
model = SPLS(n_components=1)
grid = {"penalty_x": [0.1, 0.3, 0.5, 0.9], "penalty_y": [0.1, 0.3, 0.5, 0.9]}
spls = GridSearchCV(model, grid)
plsca = PLSCanonical(n_components=1)

# Fitting SPLS model and retrieving weights
spls = spls.fit(X_train, Y_train)
spls_x_weights = spls.best_estimator_.x_weights_.T
spls_y_weights = spls.best_estimator_.y_weights_.T

# Fitting PLSCanonical model and retrieving weights
plsca.fit(X_train, Y_train)
plsca_x_weights = plsca.x_weights_.T
plsca_y_weights = plsca.y_weights_.T

# Plotting the weights of both models to visually compare their performance
# This will highlight how SPLS can effectively identify and prioritize the true signal
# over noise, demonstrating its advantage in scenarios with sparse true signals.
bar_width = 0.35  # Width of the bars

plt.figure(figsize=(12, 8))

# Plot for x-weights
plt.subplot(121)
bar1 = np.arange(len(spls_x_weights[0]))  # Positions for SPLS bars
bar2 = [x + bar_width for x in bar1]  # Positions for PLSCanonical bars

plt.bar(bar1, spls_x_weights[0], width=bar_width, label="SPLS")
plt.bar(bar2, plsca_x_weights[0], width=bar_width, label="PLSCanonical")
plt.xlabel("Variables")
plt.ylabel("Weights")
plt.title("x-weights (1st canonical vector)")
plt.xticks(
    [r + bar_width / 2 for r in range(len(spls_x_weights[0]))],
    ["Var1", "Var2", "Var3", "Var4", "Var5", "Var6"],
)
plt.legend()

# Plot for y-weights
plt.subplot(122)
bar1 = np.arange(len(spls_y_weights[0]))  # Positions for SPLS bars
bar2 = [x + bar_width for x in bar1]  # Positions for PLSCanonical bars

plt.bar(bar1, spls_y_weights[0], width=bar_width, label="SPLS")
plt.bar(bar2, plsca_y_weights[0], width=bar_width, label="PLSCanonical")
plt.xlabel("Variables")
plt.ylabel("Weights")
plt.title("y-weights (1st canonical vector)")
plt.show()
