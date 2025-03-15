"""
==================================================
Balance model complexity and cross-validated score
==================================================

This example demonstrates how to balance model complexity and cross-validated score by
finding a decent accuracy within 1 standard deviation of the best accuracy score while
minimising the number of :class:`~sklearn.decomposition.PCA` components [1]. It uses
:class:`~sklearn.model_selection.GridSearchCV with a custom refit callable to select the
optimal model.

The figure shows the trade-off between cross-validated score and the number
of PCA components. The balanced case is when `n_components=10` and `accuracy=0.88`,
which falls into the range within 1 standard deviation of the best accuracy
score.

[1] Hastie, T., Tibshirani, R.,, Friedman, J. (2001). Model Assessment and
Selection. The Elements of Statistical Learning (pp. 219-260). New York,
NY, USA: Springer New York Inc..
"""

# Authors: The scikit-learn developers
# SPDX-License-Identifier: BSD-3-Clause

# %%
# Introduction
# ------------
#
# When tuning hyperparameters, we often want to balance model complexity and
# performance. The "one-standard-error" rule is a common approach: select the simplest
# model whose performance is within one standard error of the best model's performance.
# This helps avoid overfitting by preferring simpler models when their performance is
# statistically comparable to more complex ones.

import matplotlib.pyplot as plt
import numpy as np
import polars as pl

from sklearn.datasets import load_digits
from sklearn.decomposition import PCA
from sklearn.model_selection import GridSearchCV
from sklearn.pipeline import Pipeline
from sklearn.svm import LinearSVC

# %%
# Helper functions
# ---------------
#
# We define two helper functions:
# 1. `lower_bound`: Calculates the threshold for acceptable performance
#    (best score - 1 std)
# 2. `best_low_complexity`: Selects the model with the fewest PCA components that
#    exceeds this threshold


def lower_bound(cv_results):
    """
    Calculate the lower bound within 1 standard deviation
    of the best `mean_test_scores`.

    Parameters
    ----------
    cv_results : dict of numpy(masked) ndarrays
        See attribute cv_results_ of `GridSearchCV`

    Returns
    -------
    float
        Lower bound within 1 standard deviation of the
        best `mean_test_score`.
    """
    best_score_idx = np.argmax(cv_results["mean_test_score"])

    return (
        cv_results["mean_test_score"][best_score_idx]
        - cv_results["std_test_score"][best_score_idx]
    )


def best_low_complexity(cv_results):
    """
    Balance model complexity with cross-validated score.

    Parameters
    ----------
    cv_results : dict of numpy(masked) ndarrays
        See attribute cv_results_ of `GridSearchCV`.

    Return
    ------
    int
        Index of a model that has the fewest PCA components
        while has its test score within 1 standard deviation of the best
        `mean_test_score`.
    """
    threshold = lower_bound(cv_results)
    candidate_idx = np.flatnonzero(cv_results["mean_test_score"] >= threshold)
    best_idx = candidate_idx[
        cv_results["param_reduce_dim__n_components"][candidate_idx].argmin()
    ]
    return best_idx


# %%
# Set up the pipeline and parameter grid
# -------------------------------------
#
# We create a pipeline with two steps:
# 1. Dimensionality reduction using PCA
# 2. Classification using LinearSVC
#
# We'll search over different numbers of PCA components to find the optimal complexity.

pipe = Pipeline(
    [
        ("reduce_dim", PCA(random_state=42)),
        ("classify", LinearSVC(random_state=42, C=0.01)),
    ]
)

param_grid = {"reduce_dim__n_components": [6, 8, 10, 12, 14, 16, 20, 25, 35, 45, 55]}

# %%
# Perform the search with GridSearchCV
# -----------------------------------------
#
# We use GridSearchCV with our custom `best_low_complexity` function as the refit
# parameter. This function will select the model with the fewest PCA components that
# still performs within one standard deviation of the best model.

grid = GridSearchCV(
    pipe,
    cv=10,
    n_jobs=1,  # increase this on your machine to use more physical cores
    param_grid=param_grid,
    scoring="accuracy",
    refit=best_low_complexity,
    return_train_score=True,
)

# %%
# Load the digits dataset and fit the model
# ---------------------------------------

X, y = load_digits(return_X_y=True)
grid.fit(X, y)

# %%
# Visualize the results
# -------------------
#
# We'll create a bar chart showing the test scores for different numbers of PCA
# components, along with horizontal lines indicating the best score and the
# one-standard-deviation threshold.

n_components = grid.cv_results_["param_reduce_dim__n_components"]
test_scores = grid.cv_results_["mean_test_score"]

# Create a polars DataFrame for better data manipulation and visualization
results_df = pl.DataFrame(
    {
        "n_components": n_components,
        "mean_test_score": test_scores,
        "std_test_score": grid.cv_results_["std_test_score"],
        "mean_train_score": grid.cv_results_["mean_train_score"],
        "std_train_score": grid.cv_results_["std_train_score"],
        "mean_fit_time": grid.cv_results_["mean_fit_time"],
        "rank_test_score": grid.cv_results_["rank_test_score"],
    }
)

# Sort by number of components
results_df = results_df.sort("n_components")

# Calculate the lower bound threshold
lower = lower_bound(grid.cv_results_)

# Get the best model information
best_index_ = grid.best_index_
best_components = n_components[best_index_]
best_score = grid.cv_results_["mean_test_score"][best_index_]

# Add a column to mark the selected model
results_df = results_df.with_columns(
    pl.when(pl.col("n_components") == best_components)
    .then(pl.lit("Selected"))
    .otherwise(pl.lit("Regular"))
    .alias("model_type")
)

# Get the number of CV splits from the results
n_splits = sum(
    1
    for key in grid.cv_results_.keys()
    if key.startswith("split") and key.endswith("test_score")
)

# Extract individual scores for each split
test_scores = np.array(
    [
        [grid.cv_results_[f"split{i}_test_score"][j] for i in range(n_splits)]
        for j in range(len(n_components))
    ]
)
train_scores = np.array(
    [
        [grid.cv_results_[f"split{i}_train_score"][j] for i in range(n_splits)]
        for j in range(len(n_components))
    ]
)

# Calculate mean and std of test scores
mean_test_scores = np.mean(test_scores, axis=1)
std_test_scores = np.std(test_scores, axis=1)

# Find best score and threshold
best_mean_score = np.max(mean_test_scores)
threshold = best_mean_score - std_test_scores[np.argmax(mean_test_scores)]

# Create figure with two subplots
fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(15, 6))

# Plot 1: Individual test scores with mean line and thresholds
for i, comp in enumerate(n_components):
    # Plot individual points with lower opacity
    ax1.scatter(
        [comp] * n_splits,
        test_scores[i],
        alpha=0.2,  # Lower alpha for individual points
        color="lightgray",
        label="Individual scores" if i == 0 else "",
    )

# Add mean line with error bars
ax1.errorbar(
    n_components,
    mean_test_scores,
    yerr=std_test_scores,
    color="blue",
    fmt="o-",
    capsize=5,
    label="Mean test score",
    linewidth=2,
    markersize=8,
)

# Add threshold lines
ax1.axhline(
    best_mean_score,
    color="#9b59b6",  # Purple
    linestyle="--",
    label="Best score",
    linewidth=2,
)
ax1.axhline(
    threshold,
    color="#e67e22",  # Orange
    linestyle="--",
    label="Best score - 1 std",
    linewidth=2,
)

# Highlight selected model
ax1.axvline(
    best_components,
    color="#9b59b6",  # Purple
    alpha=0.2,
    linewidth=8,
    label="Selected model",
)

# Plot 2: Train vs Test scores
for i, comp in enumerate(n_components):
    # Plot individual points
    ax2.scatter(
        [comp] * n_splits,
        train_scores[i],
        alpha=0.5,
        color="#ffd700",  # Yellow for train scores
        label="Train scores" if i == 0 else "",
    )
    ax2.scatter(
        [comp] * n_splits,
        test_scores[i],
        alpha=0.5,
        color="blue",
        label="Test scores" if i == 0 else "",
    )

# Add mean lines
ax2.plot(
    n_components,
    np.mean(train_scores, axis=1),
    "-",
    color="#ffd700",  # Yellow for train scores
    linewidth=2,
    label="Mean train score",
)
ax2.plot(
    n_components,
    np.mean(test_scores, axis=1),
    "-",
    color="blue",
    linewidth=2,
    label="Mean test score",
)

# Highlight selected model in second plot
ax2.axvline(
    best_components,
    color="#9b59b6",  # Purple
    alpha=0.2,
    linewidth=8,
    label="Selected model",
)

# Set titles and labels for both subplots
for ax in (ax1, ax2):
    ax.set_xlabel("Number of PCA components", fontsize=12)
    ax.set_xticks(n_components)
    ax.set_ylim((0.7, 1.0))
    ax.grid(True, linestyle="--", alpha=0.7)
    ax.legend(loc="lower right")

ax1.set_title("Model Selection Criteria", fontsize=14)
ax1.set_ylabel("Test accuracy", fontsize=12)

ax2.set_title("Train vs. Test Scores", fontsize=14)
ax2.set_ylabel("Score", fontsize=12)

# Add a main title for the entire figure
fig.suptitle("Balance model complexity and cross-validated score", fontsize=16, y=1.05)

# %%
# Print the results
# ---------------
#
# We print information about the selected model, including its complexity and
# performance. We also show a summary table of all models using polars.

print("Best model selected by the one-standard-error rule:")
print(f"Number of PCA components: {best_components}")
print(f"Accuracy score: {best_score:.4f}")
print(f"Best possible accuracy: {np.max(test_scores):.4f}")
print(f"Accuracy threshold (best - 1 std): {lower:.4f}")

# Create a summary table with polars
summary_df = results_df.select(
    pl.col("n_components"),
    pl.col("mean_test_score").round(4).alias("test_score"),
    pl.col("std_test_score").round(4).alias("test_std"),
    pl.col("mean_train_score").round(4).alias("train_score"),
    pl.col("std_train_score").round(4).alias("train_std"),
    pl.col("mean_fit_time").round(3).alias("fit_time"),
    pl.col("rank_test_score").alias("rank"),
)

# Add a column to mark the selected model
summary_df = summary_df.with_columns(
    pl.when(pl.col("n_components") == best_components)
    .then(pl.lit("*"))
    .otherwise(pl.lit(""))
    .alias("selected")
)

print("\nModel comparison table:")
print(summary_df)

# %%
# Conclusion
# ---------
#
# The one-standard-error rule helps us select a simpler model (fewer PCA components)
# while maintaining performance statistically comparable to the best model.
# This approach can help prevent overfitting and improve model interpretability
# and efficiency.
#
# In this example, we've seen how to implement this rule using a custom refit
# callable with :class:`~sklearn.model_selectoin.GridSearchCV`.
#
# Key takeaways:
# 1. The one-standard-error rule provides a principled way to select simpler models
# 2. Custom refit callables in :class:`~sklearn.model_selection.GridSearchCV` allow for
#    flexible model selection strategies
# 3. Visualizing both train and test scores helps identify potential overfitting
#
# This approach can be applied to other model selection scenarios where balancing
# complexity and performance is important, or in cases where a use-case specific
# selection of the "best" model is desired.

# Adjust layout and display the figure
plt.tight_layout()
plt.show()
