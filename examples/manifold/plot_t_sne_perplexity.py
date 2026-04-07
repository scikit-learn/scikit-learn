"""
=============================================================================
t-SNE: The effect of various perplexity values on the shape
=============================================================================

An illustration of t-SNE on the two concentric circles and the S-curve datasets
for different perplexity values.

We observe a tendency towards clearer shapes as the perplexity value increases.

The size, the distance and the shape of clusters may vary upon initialization,
perplexity values and does not always convey a meaning.

For further details, "How to Use t-SNE Effectively"
https://distill.pub/2016/misread-tsne/ provides a good discussion of the effects
of various parameters, as well as interactive plots to explore those effects.

"""

# Authors: The scikit-learn developers
# SPDX-License-Identifier: BSD-3-Clause

# %%
from time import time

import matplotlib.pyplot as plt
import numpy as np
from matplotlib.ticker import NullFormatter

from sklearn import datasets, manifold

n_samples = 150
perplexities = [5, 30, 50, 100]


# --- Helper ---
def plot_scatter_perplexity(ax, x, y, c, perplexity=None):
    ax.scatter(x, y, c=c)
    ax.xaxis.set_major_formatter(NullFormatter())
    ax.yaxis.set_major_formatter(NullFormatter())
    ax.axis("tight")
    if perplexity:
        ax.set_title("Perplexity=%d" % perplexity, size=14)


# --- Datasets: (label, X, colors, original_x_col, original_y_col, tsne_kwargs) ---
X_circ, y_circ = datasets.make_circles(
    n_samples=n_samples, factor=0.5, noise=0.05, random_state=0
)
X_scurve, c_scurve = datasets.make_s_curve(n_samples, random_state=0)
_t = np.linspace(0, 1, int(np.sqrt(n_samples)))
_xx, _ = np.meshgrid(_t, _t)
X_grid = np.column_stack([_xx.ravel(), _.ravel()])

dataset_configs = [
    ("Circles", X_circ, np.where(y_circ == 0, "r", "g"), X_circ[:, 0], X_circ[:, 1]),
    ("S-curve", X_scurve, c_scurve, X_scurve[:, 0], X_scurve[:, 2]),
    ("Uniform grid", X_grid, _xx.ravel(), X_grid[:, 0], X_grid[:, 1]),
]

# --- Fit & plot ---
fig, subplots = plt.subplots(3, 5, figsize=(15, 8))
fit_times = {}

tsne = manifold.TSNE(n_components=2, max_iter=300)

for row, (name, X, c, x_orig, y_orig) in enumerate(dataset_configs):
    plot_scatter_perplexity(subplots[row][0], x_orig, y_orig, c)
    fit_times[name] = []

    for col, perplexity in enumerate(perplexities, start=1):
        tsne.set_params(perplexity=perplexity)
        t0 = time()
        Y = tsne.fit_transform(X)
        elapsed = time() - t0
        fit_times[name].append(elapsed)
        print(f"{name}, perplexity={perplexity} in {elapsed:.2g} sec")
        plot_scatter_perplexity(
            subplots[row][col], Y[:, 0], Y[:, 1], c, perplexity=perplexity
        )

plt.tight_layout()
plt.show()

# %%
# For the circles dataset, t-SNE at higher perplexities recovers the two-circle
# topology, though the relative size and spacing differ from the original. The
# S-curve shape is not visually recovered at any perplexity. For the uniform
# grid, higher perplexities partially preserve the color gradient and local
# neighborhood structure, but the global geometry is distorted, illustrating
# t-SNE's tendency to prioritize local over global relationships.

# %%
fig, ax = plt.subplots(figsize=(7, 4))
for name, times in fit_times.items():
    ax.plot(perplexities, times, marker="o", label=name)
ax.set_xlabel("Perplexity")
ax.set_ylabel("Time (s)")
ax.set_title("t-SNE fit_transform time by perplexity and dataset")
ax.legend()
plt.tight_layout()
plt.show()

# %%
# The `fit_transform` time tends to increase with perplexity. Higher perplexity
# means t-SNE considers more neighbors when computing pairwise affinities,
# making that step more expensive. In this case, we only show a single
# measurement, but repeated measurements would better display the trend.
