# Authors: The scikit-learn developers
# SPDX-License-Identifier: BSD-3-Clause

"""
============================================================
Evaluating Model Performance on Highly Imbalanced Datasets
============================================================
...

In applications such as payment fraud detection or medical diagnosis, target
classes are often severely imbalanced (e.g., 1-2% positive class). In these
scenarios, standard Receiver Operating Characteristic (ROC) curves and ROC-AUC
can paint an overly optimistic picture of model performance.

This example demonstrates why Precision-Recall (PR) curves and threshold tuning
are essential for evaluating models under severe class imbalance.

"""

# Authors: scikit-learn developers
# License: BSD 3 clause

import matplotlib.pyplot as plt
import numpy as np

from sklearn.datasets import make_classification
from sklearn.ensemble import HistGradientBoostingClassifier
from sklearn.linear_model import LogisticRegression
from sklearn.metrics import PrecisionRecallDisplay, RocCurveDisplay
from sklearn.model_selection import train_test_split

# %%
# Synthetic Fraud Dataset Generation
# ----------------------------------
# We simulate a payment fraud dataset with a 99:1 class ratio (1% fraud rate)
# using :func:`~sklearn.datasets.make_classification`.

X, y = make_classification(
    n_samples=20_000,
    n_features=20,
    n_informative=12,
    n_redundant=4,
    weights=[0.99, 0.01],
    flip_y=0,
    random_state=42,
)

X_train, X_test, y_train, y_test = train_test_split(
    X, y, test_size=0.3, stratify=y, random_state=42
)

print(f"Train positive ratio: {y_train.mean():.2%}")
print(f"Test positive ratio:  {y_test.mean():.2%}")

# %%
# Model Training
# --------------
# We train a baseline :class:`~sklearn.linear_model.LogisticRegression` and a non-linear
# :class:`~sklearn.ensemble.HistGradientBoostingClassifier`.

models = {
    "Logistic Regression": LogisticRegression(max_iter=1000, random_state=42),
    "Hist Gradient Boosting": HistGradientBoostingClassifier(random_state=42),
}

for name, model in models.items():
    model.fit(X_train, y_train)

# %%
# ROC Curves vs. Precision-Recall Curves
# --------------------------------------
# Below, we plot the ROC curves and Precision-Recall curves side-by-side.

fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(13, 5))

for name, model in models.items():
    RocCurveDisplay.from_estimator(model, X_test, y_test, ax=ax1, name=name)
    PrecisionRecallDisplay.from_estimator(model, X_test, y_test, ax=ax2, name=name)

# Add baseline chance line for Precision-Recall (prevalence)
prevalence = y_test.mean()
ax2.axhline(
    prevalence,
    color="navy",
    linestyle="--",
    label=f"Chance level (prevalence = {prevalence:.2f})",
)

ax1.set_title("ROC Curves (Can be misleading under imbalance)")
ax1.grid(True, linestyle=":", alpha=0.6)

ax2.set_title("Precision-Recall Curves (Informative under imbalance)")
ax2.grid(True, linestyle=":", alpha=0.6)
ax2.legend(loc="upper right")

plt.tight_layout()
plt.show()

# %%
# Threshold Selection for Cost Trade-Offs
# --------------------------------------

hgb_model = models["Hist Gradient Boosting"]
y_probs = hgb_model.predict_proba(X_test)[:, 1]

thresholds = np.linspace(0.01, 0.99, 100)
precisions = []
recalls = []

for t in thresholds:
    y_pred_t = y_probs >= t
    tp = np.sum((y_pred_t == 1) & (y_test == 1))
    fp = np.sum((y_pred_t == 1) & (y_test == 0))
    fn = np.sum((y_pred_t == 0) & (y_test == 1))

    precision = tp / (tp + fp) if (tp + fp) > 0 else 1.0
    recall = tp / (tp + fn) if (tp + fn) > 0 else 0.0

    precisions.append(precision)
    recalls.append(recall)

fig, ax = plt.subplots(figsize=(8, 4))
ax.plot(thresholds, precisions, label="Precision", color="crimson")
ax.plot(thresholds, recalls, label="Recall", color="forestgreen")
ax.set_xlabel("Decision Threshold")
ax.set_ylabel("Score")
ax.set_title("Precision and Recall vs. Decision Threshold (HistGradientBoosting)")
ax.grid(True, linestyle=":", alpha=0.6)
ax.legend()
plt.tight_layout()
plt.show()
