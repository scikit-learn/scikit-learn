"""
=========================================================
Gradient Boosting Out-of-Bag estimates
=========================================================

This example demonstrates how to apply
:class:`~sklearn.ensemble.GradientBoostingClassifier` and monitor the
Out-of-Bag (OOB) error estimates at each iteration.
"""

# Authors: The scikit-learn developers
# SPDX-License-Identifier: BSD-3-Clause

import matplotlib.pyplot as plt
import numpy as np

from sklearn.datasets import make_hastie_10_2
from sklearn.ensemble import GradientBoostingClassifier
from sklearn.metrics import log_loss
from sklearn.model_selection import train_test_split

# %%
# Generate Data
# -------------
X, y = make_hastie_10_2(random_state=0)
X_train, X_test, y_train, y_test = train_test_split(
    X, y, test_size=0.5, random_state=42
)

# %%
# Train Gradient Boosting
# -----------------------
# We specify ``subsample=0.5`` to enable Out-of-Bag estimates.
clf = GradientBoostingClassifier(
    n_estimators=200, learning_rate=1.0, max_depth=1, random_state=0, subsample=0.5
).fit(X_train, y_train)

# %%
# Plot OOB and Test Error
# -----------------------
# We compute the test error and the OOB error at each stage of boosting.
test_score = np.zeros((clf.n_estimators,), dtype=np.float64)

for i, y_proba in enumerate(clf.staged_predict_proba(X_test)):
    test_score[i] = log_loss(y_test, y_proba)

plt.figure(figsize=(10, 6))
plt.subplot(1, 1, 1)
plt.title("Deviance")
plt.plot(
    np.arange(clf.n_estimators) + 1,
    clf.train_score_,
    "b-",
    label="Training Set Deviance",
)
plt.plot(np.arange(clf.n_estimators) + 1, test_score, "r-", label="Test Set Deviance")
plt.plot(
    np.arange(clf.n_estimators) + 1,
    clf.oob_improvement_.cumsum(),
    "k-",
    label="OOB Deviance",
)
plt.legend(loc="upper right")
plt.xlabel("Boosting Iterations")
plt.ylabel("Deviance")
plt.show()
