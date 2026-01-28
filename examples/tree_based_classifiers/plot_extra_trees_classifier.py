"""
=========================================================
Extra Trees Classifier Pixel Importance with Digits
=========================================================

This example demonstrates the usage of :class:`~sklearn.ensemble.ExtraTreesClassifier`
to visualize the importance of pixels in the Digits dataset classification.
"""

# Authors: The scikit-learn developers
# SPDX-License-Identifier: BSD-3-Clause

import matplotlib.pyplot as plt

from sklearn.datasets import load_digits
from sklearn.ensemble import ExtraTreesClassifier

# %%
# Load Data
# ---------
data = load_digits()
X, y = data.data, data.target

# %%
# Train Extra Trees Classifier
# ----------------------------
# We fit an ExtraTreesClassifier to the data. ExtraTrees are expected to be
# faster to train than Random Forests and often perform well on high dimensional
# data.
clf = ExtraTreesClassifier(n_estimators=100, random_state=0)
clf.fit(X, y)

# %%
# Visualize Pixel Importances
# ---------------------------
# We reshape the feature importances into the original image shape (8x8)
# to visualize which parts of the image are most important for classification.
importances = clf.feature_importances_
importances = importances.reshape(data.images[0].shape)

plt.figure(figsize=(8, 6))
plt.imshow(importances, cmap=plt.cm.hot, interpolation="nearest")
plt.title(
    f"Pixel importances using ExtraTreesClassifier\nAccuracy: {clf.score(X, y):.2f}"
)
plt.colorbar(label="Importance")
plt.axis("off")
plt.show()
