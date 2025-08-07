"""
===========================================================
Hierarchical clustering: structured vs unstructured Ward
===========================================================

This example demonstrates hierarchical clustering with and without
connectivity constraints. It compares standard Ward clustering with a
structured variant that enforces k-Nearest Neighbors connectivity.

Without connectivity constraints, the clustering is based purely on distance,
while with constraints, the clustering respects local structure.

In the limit of a small number of clusters, hierarchical clustering methods
tend to create a few very large clusters while leaving others almost empty.
This effect is particularly noticeable when using connectivity constraints,
as they enforce locality in cluster formation, reducing the tendency of clusters
to grow arbitrarily large. (See the discussion in Hierarchical clustering:
structured vs unstructured Ward).
"""

# Authors: The scikit-learn developers
# SPDX-License-Identifier: BSD-3-Clause

# %%
# Generate data
# -------------
# Generate the Swiss Roll dataset.
import time

from sklearn.cluster import AgglomerativeClustering
from sklearn.datasets import make_swiss_roll

n_samples = 1500
noise = 0.05
X, _ = make_swiss_roll(n_samples, noise=noise)
X[:, 1] *= 0.5  # Make the roll thinner

# Compute clustering without connectivity constraints
# ---------------------------------------------------
print("Compute unstructured hierarchical clustering...")
st = time.time()
ward_unstructured = AgglomerativeClustering(n_clusters=6, linkage="ward").fit(X)
elapsed_time_unstructured = time.time() - st
label_unstructured = ward_unstructured.labels_
print(f"Elapsed time: {elapsed_time_unstructured:.2f}s")
print(f"Number of points: {label_unstructured.size}")

# %%
# Plot unstructured clustering results
import matplotlib.pyplot as plt
import numpy as np

fig1 = plt.figure()
ax1 = fig1.add_subplot(111, projection="3d", elev=7, azim=-80)
ax1.set_position([0, 0, 0.95, 1])
for l in np.unique(label_unstructured):
    ax1.scatter(
        X[label_unstructured == l, 0],
        X[label_unstructured == l, 1],
        X[label_unstructured == l, 2],
        color=plt.cm.jet(float(l) / np.max(label_unstructured + 1)),
        s=20,
        edgecolor="k",
    )
fig1.suptitle(
    f"Without connectivity constraints (time {elapsed_time_unstructured:.2f}s)"
)

# %%
# Compute connectivity graph
# --------------------------
from sklearn.neighbors import kneighbors_graph

connectivity = kneighbors_graph(X, n_neighbors=10, include_self=False)

# %%
# Compute clustering with connectivity constraints
# ------------------------------------------------
print("Compute structured hierarchical clustering...")
st = time.time()
ward_structured = AgglomerativeClustering(
    n_clusters=6, connectivity=connectivity, linkage="ward"
).fit(X)
elapsed_time_structured = time.time() - st
label_structured = ward_structured.labels_
print(f"Elapsed time: {elapsed_time_structured:.2f}s")
print(f"Number of points: {label_structured.size}")

# %%
# Plot structured clustering results
fig2 = plt.figure()
ax2 = fig2.add_subplot(111, projection="3d", elev=7, azim=-80)
ax2.set_position([0, 0, 0.95, 1])
for l in np.unique(label_structured):
    ax2.scatter(
        X[label_structured == l, 0],
        X[label_structured == l, 1],
        X[label_structured == l, 2],
        color=plt.cm.jet(float(l) / np.max(label_structured + 1)),
        s=20,
        edgecolor="k",
    )
fig2.suptitle(f"With connectivity constraints (time {elapsed_time_structured:.2f}s)")

plt.show()
