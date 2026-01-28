"""
============================================================
Histogram-Based Gradient Boosting vs Gradient Boosting
============================================================

This example compares the training time of
:class:`~sklearn.ensemble.HistGradientBoostingClassifier` and
:class:`~sklearn.ensemble.GradientBoostingClassifier` on a larger synthetic dataset.
"""

# Authors: The scikit-learn developers
# SPDX-License-Identifier: BSD-3-Clause

import time
import matplotlib.pyplot as plt
from sklearn.datasets import make_classification
from sklearn.ensemble import GradientBoostingClassifier, HistGradientBoostingClassifier
from sklearn.model_selection import train_test_split

# %%
# Generate Large Data
# -------------------
# We generate 50,000 samples to make the speed difference noticeable.
X, y = make_classification(n_samples=50000, n_features=20, random_state=42)
X_train, X_test, y_train, y_test = train_test_split(X, y, test_size=0.2, random_state=42)

# %%
# Train Gradient Boosting
# -----------------------
print("Training GradientBoostingClassifier...")
start_time = time.time()
gb = GradientBoostingClassifier(n_estimators=50, random_state=42)
gb.fit(X_train, y_train)
gb_time = time.time() - start_time
print(f"GradientBoostingClassifier time: {gb_time:.2f}s")

# %%
# Train Hist Gradient Boosting
# ----------------------------
print("Training HistGradientBoostingClassifier...")
start_time = time.time()
hgb = HistGradientBoostingClassifier(max_iter=50, random_state=42)
hgb.fit(X_train, y_train)
hgb_time = time.time() - start_time
print(f"HistGradientBoostingClassifier time: {hgb_time:.2f}s")

# %%
# Compare Results
# ---------------
gb_score = gb.score(X_test, y_test)
hgb_score = hgb.score(X_test, y_test)

print(f"GB Accuracy: {gb_score:.4f}")
print(f"HGB Accuracy: {hgb_score:.4f}")

fig, ax = plt.subplots()
times = [gb_time, hgb_time]
names = ['GradientBoosting', 'HistGradientBoosting']
ax.bar(names, times, color=['blue', 'orange'])
ax.set_ylabel("Training Time (s)")
ax.set_title("Training Time Comparison")
plt.show()
