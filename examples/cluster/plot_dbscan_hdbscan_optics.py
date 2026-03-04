# -*- coding: utf-8 -*-
"""
=====================================================
Demo of DBSCAN, HDBSCAN, OPTICS clustering algorithms
=====================================================
.. currentmodule:: sklearn

DBSCAN, HDBSCAN, and OPTICS are density-based clustering algorithms,
meaning they leverage regional variations in density to identify
meaningful clusters. This demo will begin with DBSCAN and then move to
HDBSCAN and OPTICS to illustrate the gaps in DBSCAN that the latter
solve for.
"""

# Authors: The scikit-learn developers
# SPDX-License-Identifier: BSD-3-Clause

# %%
# Helper Functions
# ----------------
#
import matplotlib.pyplot as plt
import numpy as np


def plot(X, labels, probabilities=None, parameters=None, ground_truth=False, ax=None):
    if ax is None:
        _, ax = plt.subplots(figsize=(10, 4))
    labels = labels if labels is not None else np.ones(X.shape[0])
    probabilities = probabilities if probabilities is not None else np.ones(X.shape[0])
    # Remove black and use for noise instead.
    unique_labels = set(labels)

    colors = [plt.cm.Spectral(each) for each in np.linspace(0, 1, len(set(labels)))]
    # The probability of a point belonging to its labeled cluster
    # determines the size of its marker
    proba_map = {idx: probabilities[idx] for idx in range(len(labels))}
    for k, col in zip(unique_labels, colors):
        if k == -1:
            # Black used for noise.
            col = [0, 0, 0, 1]

        class_index = np.where(labels == k)[0]
        for ci in class_index:
            ax.plot(
                X[ci, 0],
                X[ci, 1],
                "x" if k == -1 else "o",
                markerfacecolor=tuple(col),
                markeredgecolor="black",
                alpha=0.5 if k == -1 else 0.1 + 0.9 * proba_map[ci],
            )
    n_clusters_ = len(set(labels)) - (1 if -1 in labels else 0)
    preamble = "True" if ground_truth else "Estimated"
    title = f"{preamble} number of clusters: {n_clusters_}"
    if parameters is not None:
        parameters_str = ", ".join(f"{k}={v}" for k, v in parameters.items())
        title += f" | {parameters_str}"
    ax.set_title(title)
    plt.tight_layout()


# %%
# Dataset
# -------
#
# Consider the following dataset with four clusters. An ideal
# clustering algorithm will distinguish each cluster without any prior
# information.

from sklearn.datasets import make_blobs

centers = [[-1, -1], [-1, 1], [3, 3], [3, -3]]
X, labels_true = make_blobs(
    n_samples=750, centers=centers, cluster_std=[0.2, 0.1, 0.6, 0.6], random_state=0
)

plot(X, labels=labels_true, ground_truth=True)
# %%
# DBSCAN
# ------
# DBSCAN identifies clusters by determining "core points", which are
# samples neighboring a minimum number of other samples within a
# certain distance. These two parameters are labeled `min_samples` and
# `eps` in Scikit-learn. Typically, they are not known ahead of time
# and require tuning. More details regarding the algorithm and its
# implementation can be found in :ref:`User Guide <DBSCAN>`.

# Using `eps=0.2` and `min_samples=10` on the example dataset, DBSCAN
# correctly identifies the two left-most clusters but fails for those
# on the right. This difference illustrates the limitations of choosing
# a global parameter for `min_samples` and `eps` given that the density
# varies by cluster. A key advantage of HDBSCAN and OPTICS over this
# algorithm is their ability to identify clusters at varying density
# thresholds.
from sklearn.cluster import DBSCAN

eps = 0.2
min_samples = 10

db = DBSCAN(eps=eps, min_samples=min_samples).fit(X)
plot(X, db.labels_, parameters={"eps": eps, "min_samples": min_samples})

# %%
# HDBSCAN
# -------
# HDBSCAN builds upon DBSCAN by determining clusters at varying
# densities. This is achieved by calculating the mutual reachability
# distance between pairs of data points and varying this distance. The
# full details of this algorithm can be found in :ref:`User Guide
# <HDBSCAN>`.

# Using the same dataset, HDBSCAN successfully identifies all 4 clusters
# despite the variation in density. Unlike DBSCAN and OPTICS, this
# algorithm can output probabilities for each label, which is
# demonstrated below by scaling each data point's transparency by its
# probability.

from sklearn.cluster import HDBSCAN

hdb = HDBSCAN(min_samples=min_samples, copy=False).fit(X)
plot(X, hdb.labels_, hdb.probabilities_)

# %%
# OPTICS
# ------
# Like HDBSCAN, OPTICS can be viewed as an improvement on DBSCAN, i.e.
# it generalizes the `eps` parameter to a range of values. Unlike
# HDBSCAN, however, the algorithm orders samples using the reachability
# distance. More details can be found in :ref:`User Guide<optics>`.
#
from sklearn.cluster import OPTICS, cluster_optics_dbscan

optics = OPTICS(min_samples=min_samples, min_cluster_size=0.1).fit(X)
plot(X, optics.labels_)

# %%
# Generalizations of DBSCAN
#
# HDBSCAN and OPTICS can be viewed as generalizing DBSCAN to a range of
# densities instead of a fixed value. This is illustrated below, where
# both algorithms yield similar clusters as DBSCAN when epsilon is
# fixed.

label_eps = cluster_optics_dbscan(
    reachability=optics.reachability_,
    core_distances=optics.core_distances_,
    ordering=optics.ordering_,
    eps=eps,
)

plot(X, label_eps)

# %%
label_eps = hdb.dbscan_clustering(eps)
plot(X, label_eps)
