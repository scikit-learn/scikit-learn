"""
====================================================================
Comparing StandardScaler and MinMaxScaler on Synthetic 2D Data
====================================================================

This example visually compares the effect of :class:`StandardScaler`
and :class:`MinMaxScaler` on a simple 2D dataset.

Both scalers are commonly used preprocessing techniques:

- **StandardScaler** standardizes each feature to mean 0 and std 1.
- **MinMaxScaler** rescales features to a fixed range (default: [0, 1]).

These transformations affect distance-based models such as KMeans,
KNN, and SVM. Visualizing the difference is helpful for understanding
how scaling changes the geometry of the data.
"""

import numpy as np
import matplotlib.pyplot as plt
from sklearn.preprocessing import StandardScaler, MinMaxScaler

# -----------------------------
# Generate synthetic dataset
# -----------------------------
rng = np.random.RandomState(42)
X = np.vstack([
    rng.normal(loc=[2, 8], scale=[1.0, 2.0], size=(200, 2)),
    rng.normal(loc=[8, 3], scale=[2.0, 1.0], size=(200, 2)),
])

# -----------------------------
# Apply the scalers
# -----------------------------
std_scaler = StandardScaler()
mm_scaler = MinMaxScaler()

X_std = std_scaler.fit_transform(X)
X_mm = mm_scaler.fit_transform(X)

# -----------------------------
# Plotting function
# -----------------------------
def plot(ax, data, title):
    ax.scatter(data[:, 0], data[:, 1], s=15, alpha=0.6)
    ax.set_title(title)
    ax.set_xlabel("Feature 1")
    ax.set_ylabel("Feature 2")


# -----------------------------
# Create plots
# -----------------------------
fig, axes = plt.subplots(1, 3, figsize=(15, 4))

plot(axes[0], X, "Original Data")
plot(axes[1], X_std, "After StandardScaler")
plot(axes[2], X_mm, "After MinMaxScaler")

fig.suptitle("Effect of Feature Scaling")
plt.tight_layout()
plt.show()

# -----------------------------
# Summary
# -----------------------------
#
# StandardScaler centers and scales features independently, preserving
# relative distances but changing absolute coordinates.
#
# MinMaxScaler maps each feature to a fixed range, which can distort
# distances if outliers are present.
#
# Choosing the right scaler depends on the downstream model:
#
# - Distance-based models (KMeans, KNN) benefit from StandardScaler.
# - Neural networks often benefit from MinMaxScaler.
#

