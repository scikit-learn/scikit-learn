"""
======================================================================
Benchmark: Comparing Isomap Solvers - Execution Time vs. Representation
======================================================================

This benchmark demonstrates how different eigenvalue solvers in Isomap 
can affect execution time and embedding quality.

Description:
------------
We use a subset of handwritten digits (`load_digits` with 6 classes). 
Each data point is projected into a lower-dimensional space (2D) using 
two different solvers (`auto` and `randomized`).

What you can observe:
----------------------
- The `auto` solver provides a reference solution.
- The `randomized` solver is tested for comparison in terms of 
  representation quality and execution time.

Further exploration:
---------------------
You can modify the number of neighbors (`n_neighbors`) or experiment with 
other Isomap solvers.
"""

import time
import numpy as np
import matplotlib.pyplot as plt
from matplotlib import offsetbox
from sklearn.datasets import load_digits
from sklearn.preprocessing import MinMaxScaler
from sklearn.manifold import Isomap

# 1- Data Loading
# ---------------
digits = load_digits(n_class=6)
X, y = digits.data, digits.target
n_neighbors = 30  # Number of neighbors for Isomap

# 2- Visualization Function
# -------------------------
def plot_embedding(ax, X, title):
    """Displays projected points with image annotations."""
    X = MinMaxScaler().fit_transform(X)

    for digit in digits.target_names:
        ax.scatter(
            *X[y == digit].T,
            marker=f"${digit}$",
            s=60,
            color=plt.cm.Dark2(digit),
            alpha=0.425,
            zorder=2,
        )

    # Add digit images in the projected space
    shown_images = np.array([[1.0, 1.0]])
    for i in range(X.shape[0]):
        dist = np.sum((X[i] - shown_images) ** 2, 1)
        if np.min(dist) < 4e-3:
            continue
        shown_images = np.concatenate([shown_images, [X[i]]], axis=0)
        imagebox = offsetbox.AnnotationBbox(
            offsetbox.OffsetImage(digits.images[i], cmap=plt.cm.gray_r), X[i]
        )
        imagebox.set(zorder=1)
        ax.add_artist(imagebox)

    ax.set_title(title)
    ax.axis("off")

# 3- Define Embeddings and Benchmark
# ----------------------------------
embeddings = {
    "Isomap (auto solver)": Isomap(n_neighbors=n_neighbors, n_components=2, eigen_solver='auto'),
    "Isomap (randomized solver)": Isomap(n_neighbors=n_neighbors, n_components=2, eigen_solver='randomized_value'),
}

projections, timing = {}, {}

# Compute embeddings
for name, transformer in embeddings.items():
    print(f"Computing {name}...")
    start_time = time.time()
    projections[name] = transformer.fit_transform(X, y)
    timing[name] = time.time() - start_time

# 4- Display Results
# ------------------
fig, axes = plt.subplots(1, 2, figsize=(14, 6))

for ax, (name, proj) in zip(axes, projections.items()):
    title = f"{name} (time: {timing[name]:.3f}s)"
    plot_embedding(ax, proj, title)

plt.tight_layout()
plt.show()
