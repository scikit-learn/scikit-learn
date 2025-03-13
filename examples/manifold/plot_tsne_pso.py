"""Comparison of t-SNE and t-SNE with Particle Swarm Optimization (PSO)."""

# Authors: The scikit-learn developers
# SPDX-License-Identifier: BSD-3-Clause

import time

import matplotlib.pyplot as plt
import numpy as np
from matplotlib.ticker import NullFormatter

from sklearn import datasets
from sklearn.cluster import KMeans

# Next line to silence pyflakes. This import is needed.
from sklearn.manifold import TSNE, TSNEPSO
from sklearn.metrics import (
    adjusted_rand_score,
    normalized_mutual_info_score,
    silhouette_score,
)
from sklearn.preprocessing import StandardScaler

print(__doc__)

# ----------------------------------------------------------------------
# Load datasets

# Load Digits dataset
digits = datasets.load_digits(n_class=6)
X_digits = digits.data
y_digits = digits.target

# Load Iris dataset
iris = datasets.load_iris()
X_iris = iris.data
y_iris = iris.target

# ----------------------------------------------------------------------
# Scale data

X_digits = StandardScaler().fit_transform(X_digits)
X_iris = StandardScaler().fit_transform(X_iris)

# ----------------------------------------------------------------------
# Define evaluation metrics for comparing embeddings


def evaluate_embedding(X_embedded, y_true, n_clusters):
    """Evaluate embedding quality using clustering metrics."""
    # Apply KMeans to the embedding
    kmeans = KMeans(n_clusters=n_clusters, n_init=10, random_state=0)
    y_pred = kmeans.fit_predict(X_embedded)

    # Calculate metrics
    silhouette = silhouette_score(X_embedded, y_pred)
    nmi = normalized_mutual_info_score(y_true, y_pred)
    ari = adjusted_rand_score(y_true, y_pred)

    return {"silhouette": silhouette, "nmi": nmi, "ari": ari}


# ----------------------------------------------------------------------
# Perform dimensionality reduction on digits

# Common parameters
n_components = 2
perplexity = 30.0
random_state = 0

# Run t-SNE
print("Computing t-SNE embedding for digits")
t0 = time.time()
tsne = TSNE(
    n_components=n_components,
    perplexity=perplexity,
    n_iter=1000,
    random_state=random_state,
)
X_digits_tsne = tsne.fit_transform(X_digits)
print(f"t-SNE done in {time.time() - t0:.3f} seconds")

# Run t-SNE-PSO
print("Computing t-SNE-PSO embedding for digits")
t0 = time.time()
tsne_pso = TSNEPSO(
    n_components=n_components,
    perplexity=perplexity,
    n_iter=1000,
    n_particles=10,
    use_hybrid=True,
    random_state=random_state,
)
X_digits_tsne_pso = tsne_pso.fit_transform(X_digits)
print(f"t-SNE-PSO done in {time.time() - t0:.3f} seconds")
print(f"t-SNE-PSO final KL divergence: {tsne_pso.kl_divergence_:.4f}")

# Evaluate embeddings
n_classes_digits = len(np.unique(y_digits))
digits_tsne_metrics = evaluate_embedding(X_digits_tsne, y_digits, n_classes_digits)
digits_tsne_pso_metrics = evaluate_embedding(
    X_digits_tsne_pso, y_digits, n_classes_digits
)

print("\nDigits dataset metrics:")
print(
    f"t-SNE - Silhouette: {digits_tsne_metrics['silhouette']:.4f}, "
    f"NMI: {digits_tsne_metrics['nmi']:.4f}, "
    f"ARI: {digits_tsne_metrics['ari']:.4f}"
)
print(
    f"t-SNE-PSO - Silhouette: {digits_tsne_pso_metrics['silhouette']:.4f}, "
    f"NMI: {digits_tsne_pso_metrics['nmi']:.4f}, "
    f"ARI: {digits_tsne_pso_metrics['ari']:.4f}"
)

# ----------------------------------------------------------------------
# Perform dimensionality reduction on iris

# Run t-SNE
print("\nComputing t-SNE embedding for iris")
t0 = time.time()
tsne = TSNE(
    n_components=n_components,
    perplexity=perplexity,
    n_iter=1000,
    random_state=random_state,
)
X_iris_tsne = tsne.fit_transform(X_iris)
print(f"t-SNE done in {time.time() - t0:.3f} seconds")

# Run t-SNE-PSO
print("Computing t-SNE-PSO embedding for iris")
t0 = time.time()
tsne_pso = TSNEPSO(
    n_components=n_components,
    perplexity=perplexity,
    n_iter=1000,
    n_particles=10,
    use_hybrid=True,
    random_state=random_state,
)
X_iris_tsne_pso = tsne_pso.fit_transform(X_iris)
print(f"t-SNE-PSO done in {time.time() - t0:.3f} seconds")
print(f"t-SNE-PSO final KL divergence: {tsne_pso.kl_divergence_:.4f}")

# Evaluate embeddings
n_classes_iris = len(np.unique(y_iris))
iris_tsne_metrics = evaluate_embedding(X_iris_tsne, y_iris, n_classes_iris)
iris_tsne_pso_metrics = evaluate_embedding(X_iris_tsne_pso, y_iris, n_classes_iris)

print("\nIris dataset metrics:")
print(
    f"t-SNE - Silhouette: {iris_tsne_metrics['silhouette']:.4f}, "
    f"NMI: {iris_tsne_metrics['nmi']:.4f}, "
    f"ARI: {iris_tsne_metrics['ari']:.4f}"
)
print(
    f"t-SNE-PSO - Silhouette: {iris_tsne_pso_metrics['silhouette']:.4f}, "
    f"NMI: {iris_tsne_pso_metrics['nmi']:.4f}, "
    f"ARI: {iris_tsne_pso_metrics['ari']:.4f}"
)

# ----------------------------------------------------------------------
# Visualize the results


def plot_embedding(X, y, title, ax):
    """Plot an embedding X with the class labels y."""
    x_min, x_max = np.min(X, 0), np.max(X, 0)
    X = (X - x_min) / (x_max - x_min)

    for i in range(X.shape[0]):
        ax.text(
            X[i, 0],
            X[i, 1],
            str(y[i]),
            color=plt.cm.tab10(y[i] / 10.0),
            fontdict={"weight": "bold", "size": 9},
        )

    ax.set_title(title)
    ax.set_xticks([])
    ax.set_yticks([])
    return ax


# Create the figure
fig, axes = plt.subplots(2, 2, figsize=(15, 10))

# Plot digits embeddings
plot_embedding(
    X_digits_tsne, y_digits, "t-SNE embedding of the digits dataset", axes[0, 0]
)
plot_embedding(
    X_digits_tsne_pso, y_digits, "t-SNE-PSO embedding of the digits dataset", axes[0, 1]
)

# Plot iris embeddings
plot_embedding(X_iris_tsne, y_iris, "t-SNE embedding of the iris dataset", axes[1, 0])
plot_embedding(
    X_iris_tsne_pso, y_iris, "t-SNE-PSO embedding of the iris dataset", axes[1, 1]
)

plt.tight_layout()
plt.show()

# ----------------------------------------------------------------------
# Scatter plots showing cluster separation

fig, axes = plt.subplots(2, 2, figsize=(15, 10))

# Add colors for visualization
colors = plt.cm.tab10(np.linspace(0, 1, 10))

# Plot digits embeddings
axes[0, 0].scatter(
    X_digits_tsne[:, 0],
    X_digits_tsne[:, 1],
    c=y_digits,
    cmap=plt.cm.tab10,
    edgecolor="k",
    s=50,
)
axes[0, 0].set_title("t-SNE embedding of the digits dataset")
axes[0, 0].xaxis.set_major_formatter(NullFormatter())
axes[0, 0].yaxis.set_major_formatter(NullFormatter())

axes[0, 1].scatter(
    X_digits_tsne_pso[:, 0],
    X_digits_tsne_pso[:, 1],
    c=y_digits,
    cmap=plt.cm.tab10,
    edgecolor="k",
    s=50,
)
axes[0, 1].set_title("t-SNE-PSO embedding of the digits dataset")
axes[0, 1].xaxis.set_major_formatter(NullFormatter())
axes[0, 1].yaxis.set_major_formatter(NullFormatter())

# Plot iris embeddings
axes[1, 0].scatter(
    X_iris_tsne[:, 0],
    X_iris_tsne[:, 1],
    c=y_iris,
    cmap=plt.cm.tab10,
    edgecolor="k",
    s=50,
)
axes[1, 0].set_title("t-SNE embedding of the iris dataset")
axes[1, 0].xaxis.set_major_formatter(NullFormatter())
axes[1, 0].yaxis.set_major_formatter(NullFormatter())

axes[1, 1].scatter(
    X_iris_tsne_pso[:, 0],
    X_iris_tsne_pso[:, 1],
    c=y_iris,
    cmap=plt.cm.tab10,
    edgecolor="k",
    s=50,
)
axes[1, 1].set_title("t-SNE-PSO embedding of the iris dataset")
axes[1, 1].xaxis.set_major_formatter(NullFormatter())
axes[1, 1].yaxis.set_major_formatter(NullFormatter())

plt.tight_layout()
plt.show()
