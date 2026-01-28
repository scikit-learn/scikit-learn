"""
=========================================================
Random Forest Classifier Feature Importance
=========================================================

This example demonstrates the usage of :class:`~sklearn.ensemble.RandomForestClassifier`
to evaluate the importance of features in an artificial classification task.
"""

# Authors: The scikit-learn developers
# SPDX-License-Identifier: BSD-3-Clause

import matplotlib.pyplot as plt
import numpy as np

from sklearn.datasets import make_classification
from sklearn.ensemble import RandomForestClassifier
from sklearn.model_selection import train_test_split

# %%
# Generate Data
# -------------
# We create a synthetic dataset with 10 features, of which 3 are informative,
# 2 are redundant, and the rest are noise.
X, y = make_classification(
    n_samples=1000,
    n_features=10,
    n_informative=3,
    n_redundant=2,
    n_repeated=0,
    n_classes=2,
    random_state=0,
    shuffle=False,
)

# Split into train and test
X_train, X_test, y_train, y_test = train_test_split(
    X, y, stratify=y, random_state=42
)

# %%
# Train Random Forest
# -------------------
feature_names = [f"feature {i}" for i in range(X.shape[1])]
forest = RandomForestClassifier(random_state=0)
forest.fit(X_train, y_train)

# %%
# Feature Importance
# ------------------
# Feature importances are provided by the fitted attribute ``feature_importances_``
# and they are computed as the mean and standard deviation of accumulation of
# the impurity decrease within each tree.
importances = forest.feature_importances_
std = np.std([tree.feature_importances_ for tree in forest.estimators_], axis=0)

# %%
# Plot Importances
# ----------------
import pandas as pd

forest_importances = pd.Series(importances, index=feature_names)

fig, ax = plt.subplots()
forest_importances.plot.bar(yerr=std, ax=ax)
ax.set_title("Feature importances using MDI")
ax.set_ylabel("Mean decrease in impurity")
fig.tight_layout()
plt.show()

print(f"Test Set Accuracy: {forest.score(X_test, y_test):.2f}")
