"""
============================================
Evaluating a ranking model with NDCG, Kendall's tau and Spearman's rho
============================================

This example trains a simple model on a ranking-style problem and shows
how to evaluate it with three different metrics: NDCG, Kendall's tau and
Spearman's rho. It also explains when you would reach for each one.
"""

# Authors: The scikit-learn developers
# SPDX-License-Identifier: BSD-3-Clause

# %%
# Create a ranking-style dataset
# -------------------------------
# In a ranking problem, items are grouped into queries, and within each
# query the items need to be sorted by relevance. Here we simulate this by
# generating a regression dataset and splitting it into groups, where each
# group stands in for one query and the target is the true relevance score.

import numpy as np

from sklearn.datasets import make_regression

n_queries = 40
n_items_per_query = 10

X, y = make_regression(
    n_samples=n_queries * n_items_per_query,
    n_features=5,
    noise=20,
    random_state=0,
)

# NDCG expects relevance scores to be non-negative, just like in a real
# ranking dataset, so we shift the target above zero.
y = y - y.min()

X_groups = X.reshape(n_queries, n_items_per_query, -1)
y_groups = y.reshape(n_queries, n_items_per_query)

# %%
# Split into train and test queries
# ------------------------------------
# We hold out entire queries for testing, rather than individual items, so
# that we evaluate the model on queries it has never seen before.

from sklearn.model_selection import train_test_split

train_idx, test_idx = train_test_split(
    np.arange(n_queries), test_size=0.3, random_state=0
)

X_train = X_groups[train_idx].reshape(-1, X.shape[1])
y_train = y_groups[train_idx].reshape(-1)

X_test = X_groups[test_idx]
y_true_groups = y_groups[test_idx]

# %%
# Train a model
# --------------
# We fit a regressor on the training queries to predict a relevance score
# for each item. In a real ranking problem, items would then be sorted by
# this predicted score. We evaluate on the held-out test queries.

from sklearn.ensemble import GradientBoostingRegressor

model = GradientBoostingRegressor(random_state=0)
model.fit(X_train, y_train)

n_test_queries = len(test_idx)
X_test_flat = X_test.reshape(-1, X.shape[1])
y_pred_groups = model.predict(X_test_flat).reshape(n_test_queries, n_items_per_query)

# %%
# Evaluate with NDCG
# -------------------
# NDCG (Normalized Discounted Cumulative Gain) checks whether the top items
# in the predicted ranking are actually the most relevant ones. Mistakes
# near the top of the list hurt the score much more than mistakes near the
# bottom, which makes it a good fit for problems like search results, where
# users mostly look at the first few items.

from sklearn.metrics import ndcg_score

ndcg = ndcg_score(y_true_groups, y_pred_groups)
print(f"NDCG score: {ndcg:.3f}")

# %%
# Evaluate with Kendall's tau and Spearman's rho
# ------------------------------------------------
# Kendall's tau and Spearman's rho both measure how closely two full
# rankings agree, treating every position in the list equally. A score of
# 1 means the rankings match perfectly, 0 means no relationship, and -1
# means they are completely reversed.
#
# scipy returns an object with both a statistic and a p-value, so we only
# keep the statistic here.

from scipy.stats import kendalltau, spearmanr

kendall_scores = []
spearman_scores = []

for true_row, pred_row in zip(y_true_groups, y_pred_groups):
    kendall_scores.append(kendalltau(true_row, pred_row).statistic)
    spearman_scores.append(spearmanr(true_row, pred_row).statistic)

print(f"Average Kendall's tau: {np.mean(kendall_scores):.3f}")
print(f"Average Spearman's rho: {np.mean(spearman_scores):.3f}")

# %%
# Visualize the results
# -----------------------
# The left plot shows true relevance against predicted relevance for one
# test query. Points close to a straight diagonal line mean the model
# ranked that query well. The right plot compares the three metric scores
# side by side.

import matplotlib.pyplot as plt

fig, axes = plt.subplots(1, 2, figsize=(10, 4))

axes[0].scatter(y_true_groups[0], y_pred_groups[0])
axes[0].set_xlabel("True relevance")
axes[0].set_ylabel("Predicted relevance")
axes[0].set_title("True vs predicted relevance\n(first test query)")

metric_names = ["NDCG", "Kendall's tau", "Spearman's rho"]
metric_values = [ndcg, np.mean(kendall_scores), np.mean(spearman_scores)]
axes[1].bar(metric_names, metric_values)
axes[1].set_ylim(0, 1)
axes[1].set_title("Metric scores")

fig.tight_layout()
plt.show()

# %%
# When to use which metric
# --------------------------
# Use NDCG when only the top of the list matters, such as a search engine
# results page, where users rarely scroll past the first few results.
#
# Use Kendall's tau or Spearman's rho when the full ranking matters, not
# just the top, such as ranking every candidate in a hiring pipeline, where
# a mistake anywhere in the list is equally costly.
