# -*- coding: utf-8 -*-
"""
=======================================================
Demo of DBSCAN, HDBSCAN, OPTICS clustering algorithms
=======================================================
.. currentmodule:: sklearn

In this demo we will take a look at DBSCAN, HDBSCAN, and OPTICs clustering
algorithms. We will run each algorithm on multiclass datasets of varying
densities, and note the performance changes as we adjust important
hyperparameters. We will also show how OPTICS and HDBSCAN
can be viewed as generalizations of DBSCAN, and efficiently extract DBSCAN
clusterings from the results of running these algorithms.

We start by defining helper functions to visualize a dataset and the resulting
cluster labels after running a clustering algorithm, if applicable. We will
use this function throughout this document.
"""

# Authors: The scikit-learn developers
# SPDX-License-Identifier: BSD-3-Clause

# %%
import matplotlib.pyplot as plt
import numpy as np


def cluster_colours(labels):
    return [plt.cm.Spectral(each) for each in np.linspace(0, 1, len(set(labels)))]


def plot(X, labels, probabilities=None, parameters=None, ground_truth=False, ax=None):
    if ax is None:
        _, ax = plt.subplots(figsize=(10, 4))
    labels = labels if labels is not None else np.ones(X.shape[0])
    probabilities = probabilities if probabilities is not None else np.ones(X.shape[0])
    # Black removed and is used for noise instead.
    unique_labels = set(labels)

    colors = cluster_colours(labels)
    # The probability of a point belonging to its labeled cluster determines
    # the size of its marker
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
                markersize=4 if k == -1 else 1 + 5 * proba_map[ci],
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
# DBSCAN
# -------
# DBSCAN stands for Density-Based Spatial Clustering of Applications with Noise.
# It is a clustering algorithm that takes in a dataset of points and a measure
# of distance between them.
#
# To gain some intuition about how DBSCAN works, we
# need to discuss the concept of core points. We will also discuss neighborhoods.
# Generally speaking, when we talk about the neighborhood of a point, we are
# referring to a collection of points in the dataset that are within a certain
# distance from it. The size of the neighborhood will depend on the context. A
# small neighborhood will contain only points that are very close to the point
# in question, while a larger neighborhood will contain points that are further
# away.
#
# Core points are points in the dataset with 'enough' points near them. In
# other words, a core point is a point whose neighborhood has a minimum number
# of points in it. Clusters in DBSCAN are assigned so that any point in a
# neighborhood of a core point is considered to be in the same cluster as that
# core point.
#
# To formalize these notions, DBSCAN has two parameters that fully
# define what makes a sample a core point. These are `eps` and `min_samples`.
#
# `eps` is a positive distance. For a sample to be a core point, there must be
# `min_samples` samples in the dataset that are within `eps` distance from that
# sample. In other words, the neighborhood of a core point contains all points
# that are within `eps` distance from it, and this collection must contain at
# least `min_samples` points.
#
# Let us now try running some examples. First we will need to generate a
# simple ground-truth dataset. We use `make_blobs` to generate a dataset
# from three two-dimensional isotropic Gaussian distributions. For ease of
# presentation, our distance metric throughout this document will be Euclidian
# distance of points in the two-dimensional plane (all clustering algorithms
# discussed here can be given an optional metric to be used in place of this one).

from sklearn.datasets import make_blobs

centers = [[-1, -1], [-1, 1], [1, 0]]
X, labels_true = make_blobs(
    n_samples=750, centers=centers, cluster_std=[0.2, 0.35, 0.5], random_state=0
)

plot(X, labels=labels_true, ground_truth=True)
# %%
# Now let us see how DBSCAN performs on this dataset. We will use
# `eps=0.2` and `min_samples=10`. Very shortly we will discuss how one
# chooses parameter values according to the dataset they work with.
from sklearn.cluster import DBSCAN

db = DBSCAN(eps=0.2, min_samples=10).fit(X)
plot(X, db.labels_, parameters={"eps": 0.2, "min_samples": 10})

# %%
# We see that with these settings, DBSCAN correctly computed the number of
# clusters, and labelled most points accurately. Note that many points in the
# two sparser clusters were labelled as noise, because they were not within
# reachable from a core point (the euclidian distance from a core point was at
# least `eps`).
#
# Tuning DBSCAN
# ++++++++++++++
# How do we choose what parameter values to use? While DBSCAN provides
# default values for `eps` and `min_samples` they usually will not provide good
# results, as these sensitive to the shape of the dataset you are working with.
# Larger values of `min_samples` yield more robustness to noise, but
# with a risk of grouping small clusters together. Smaller values will lead
# to more clusters and less noise. The same is true for the `eps` parameter.
# Typically `eps` is the more sensitive parameter, so one may tune `eps` after
# finding a good value for `min_samples`.
#
# We illustrate this by running DBSCAN for `eps` values 0.1, 0.3, 0.4, while
# leaving `min_samples` fixed to 10. For consistency, we will continue to use
# these values for `eps` and `min_samples` when applicable.
eps_values = [0.1, 0.3, 0.4]

fig, axes = plt.subplots(len(eps_values), 1, figsize=(10, 12))

for idx, eps in enumerate(eps_values):
    db = DBSCAN(eps=eps, min_samples=10).fit(X)
    plot(X, db.labels_, parameters={"eps": eps, "min_samples": 10}, ax=axes[idx])

# %%
# DBSCAN produces the best results when we use the parameters as in the second
# plot. However, finding good values of `eps` and `min_samples` is not always
# as easy as this. It may require specialized knowledge of the dataset
# one is working with. While standardizing the dataset with
# `StandardScaler` may help with this problem, great care must be taken for
# choosing the correct value of `eps`.
#
# To illustrate these difficulties, we generate another dataset in the same
# way as before, but with more variability in the clusters.

centers = [[-0.85, -0.85], [-0.85, 0.85], [3, 3], [3, -3]]
X, labels_true = make_blobs(
    n_samples=750, centers=centers, cluster_std=[0.2, 0.35, 1.35, 1.35], random_state=0
)
plot(X, labels=labels_true, ground_truth=True)

# %%
# Now we run DBSCAN on this dataset with the same settings as before.
_, axes = plt.subplots(len(eps_values), 1, figsize=(10, 12))

for idx, eps in enumerate(eps_values):
    db = DBSCAN(eps=eps, min_samples=10).fit(X)
    plot(X, db.labels_, parameters={"eps": eps, "min_samples": 10}, ax=axes[idx])
# %%
# We can see that DBSCAN has a lot more trouble with the multiscale clusters
# of this dataset. For smaller values of `eps`, we label too many points as noise,
# especially in the sparse clusters. But taking `eps` as 0.4 results in
# classifying the two smaller dense clusters as one.
#
# The difficulty of DBSCAN with multiscale datasets can be attributed to the
# fixed value of `eps` in the definition of core points in the DBSCAN algorithm,
# resulting in fixed neighborhood sizes when classifying the core points in the
# dataset. In contrast to this approach, HDBSCAN and OPTICS generalize the `eps`
# parameter in DBSCAN to a range of values. Generally speaking they have better
# performance on multiscale clusters.
#
# Before we move on to discuss these algorithms we first discuss evaluating
# clustering performance. The dataset you have might not be the easiest to
# visualize in the way we are doing. In this case, another option to evaluate
# performance of a clustering algorithm is to use metrics.

# %%
# Measuring Clustering Performance
# -----------------------------------
# If the ground-truth labels are known, as in our current situation, we can use
# metrics that quantify the quality of the resulting clusters. Examples
# include homogeneity, completeness, V-measure, Rand-Index,
# Adjusted Rand-Index and Adjusted Mutual Information (AMI).
#
# If the ground-truth labels are not known, we only have the results of running
# the clustering algorithm. In this case we may use the Silhouette Coefficient.
# For brevity, we will be satisfied with a simple demonstration of evaluating
# these metrics with the result of our last DBSCAN run.
#
# For more information, see the
# :ref:`sphx_glr_auto_examples_cluster_plot_adjusted_for_chance_measures.py`
# example or the :ref:`clustering_evaluation` module.
from sklearn import metrics

labels = db.labels_

print(f"Homogeneity: {metrics.homogeneity_score(labels_true, labels):.3f}")
print(f"Completeness: {metrics.completeness_score(labels_true, labels):.3f}")
print(f"V-measure: {metrics.v_measure_score(labels_true, labels):.3f}")
print(f"Adjusted Rand Index: {metrics.adjusted_rand_score(labels_true, labels):.3f}")
print(
    "Adjusted Mutual Information:"
    f" {metrics.adjusted_mutual_info_score(labels_true, labels):.3f}"
)

print(f"Silhouette Coefficient: {metrics.silhouette_score(X, labels):.3f}")

# %%
# HDBSCAN
# ---------
# We now move onto HDBSCAN, or hierarchical DBSCAN. This algorithm can be
# viewed as an improvement to the ideas in DBSCAN, but takes into account
# different scales of the dataset, using the notion of core-distances and
# reachability-distances instead of core points with fixed `eps` sized
# neighborhoods. For brevity we will leave the full details of this algorithm
# for the :ref:`User Guide <HDBSCAN>`.
#
# One of the greatest advantages of HDBSCAN over DBSCAN is its out-of-the-box
# robustness. It’s especially remarkable on heterogeneous mixtures of data. For
# a quick example, let us run it on the previous dataset with its default
# settings.
from sklearn.cluster import HDBSCAN

hdb = HDBSCAN().fit(X)
plot(X, hdb.labels_, hdb.probabilities_)

# %%
# HDBSCAN is able to adapt to the multi-scale structure of the dataset without
# requiring any parameter tuning. Also notice that on fit, HDBSCAN produces a
# probabilities attribute that indicates, for each sample, the strength that it
# is a member of its assigned cluster. Our plotting function takes in this
# attribute to enlargen marker size of a sample based on its strength value.
#
# Tuning HDBSCAN
# ++++++++++++++++
# HDBSCAN still has some hyperparemeters, and a more interesting dataset will
# probably require tuning them. We focus on two of them, `min_cluster_size` and
# `min_samples`
#
# `min_cluster_size` is the minimum number of samples in a group for that group
# to be considered a cluster.
#
# This is an intuitive parameter to select for. Clusters smaller than the ones
# of this size will be left as noise. The default value is 5. This parameter is
# generally tuned to larger values as needed. Smaller values will likely to
# lead to results with fewer points labeled as noise. However values which too
# small will lead to false sub-clusters being picked up and preferred. Larger
# values tend to be more robust with respect to noisy datasets, e.g.
# high-variance clusters with significant overlap.
#
# `min_samples` corresponds to, in some sense, the number of points required in
# the neighborhood of a point, for that point to be considered a core point. It
# is similar to its counterpart in DBSCAN.
#
# A full definition of this parameter requires more context, which the reader
# can find in the :ref:`User Guide <HDBSCAN>`. For our purposes, one can
# interpret this parameter as a measure of how conservative we want our
# clustering to be. Larger values for `min_samples` increase the model’s
# robustness to noise, but risks ignoring or discarding potentially valid but
# small clusters.
#
# We recommend tuning `min_samples` after finding a good value for
# `min_cluster_size`. Alternatively one can set `min_samples=min_cluster_size`
# and simplify the hyperparemter space. For brevity, this is the option we are
# going to take in the example plots below as we show the results of running
# HDBSCAN with different hyperparameter values.

_, axes = plt.subplots(len(eps_values), 1, figsize=(10, 12))
hdbscan_param_vals = [5, 15, 25]
for idx, val in enumerate(hdbscan_param_vals):
    # default value for min_samples is min_cluster_size but
    # we set it explicitly for clarity
    hdb = HDBSCAN(min_cluster_size=val, min_samples=val).fit(X)
    plot(
        X,
        hdb.labels_,
        hdb.probabilities_,
        parameters={"min_cluster_size": val, "min_samples": val},
        ax=axes[idx],
    )


# %%
# Extracting DBSCAN clusters
# ++++++++++++++++++++++++++
# During fit, HDBSCAN builds a (singkle-linkage) tree which encodes a DBSCAN
# clustering of all points across all values of DBSCAN’s `eps` parameter. We
# can efficiently obtain these DBSCAN like clusterings efficiently without
# fully recomputing intermediate values required for an HDBSCAN fit, such as
# core-distances, mutual-reachability, and the minimum spanning tree This is
# done by specifying a `cut_distance` (equivalent to `eps`) that we want to
# cluster with. Again, we refer to the :ref:`User Guide <HDBSCAN>` for more
# details.
#
# We run HDBSCAN with `min_cluster_size` and `min_samples` set to 10 and then
# use `HDBSCAN.dbscan_clustering(eps)` for various `eps` (the ones we used
# previously when discussing DBSCAN). Compare the plot generated below with the
# one we previously generated by running DBSCAN on the same dataset, and note
# the similarity.

hdb = HDBSCAN(min_cluster_size=10, min_samples=10)
hdb.fit(X)
fig, axes = plt.subplots(len(eps_values), 1, figsize=(10, 12))
for idx, eps in enumerate(eps_values):
    params = {"threshold": eps}
    hdb_dbscan_clustering_labels = hdb.dbscan_clustering(eps)
    plot(
        X,
        hdb_dbscan_clustering_labels,
        hdb.probabilities_,
        parameters=params,
        ax=axes[idx],
    )

# %%
# OPTICS
# -------
# We now discuss OPTICS (Ordering Points To Identify the Clustering Structure).
# Like HDBSCAN, OPTICS can be viewed as an improvement on DBSCAN. It
# generalizes the `eps` parameter from DBSCAN to a range of values. To specify
# this range, the user can set the parameter called `max_eps` which OPTICS uses
# when looking at the sizes of neighborhoods of points. The default value is
# `np.inf`. Smaller values of `max_eps` will result in shorter runtimes.
#
# It also has parametmers `min_samples` and `min_cluster_size`, whose
# descriptions are similar to their counterparts in HDBSCAN, so we will not
# repeat them here. It should be noted that values for these parameters that
# work well for HDBSCAN might not work well for OPTICS, and vice versa.
#
# `min_cluster_size` can either be a positive integer or a float between 0 and
# 1. In the latter case, it is specified as a fraction of the total number of
# samples in the dataset.
from sklearn.cluster import OPTICS, cluster_optics_dbscan

# min_cluster_size is 10% of len(X)
optics = OPTICS(min_samples=10, min_cluster_size=0.1).fit(X)
plot(X, optics.labels_)

# %%
# A key difference between OPTICS and the former algorithms is that OPTICS
# assigns each sample a reachability distance and an ordering. The assigned
# reachability distance for a sample can be thought of, informally, as the
# distance to the closest core point. The ordering in which OPTICS labels the
# points is such that points that are close together in distance (in the sense
# of the metric used by OPTICS, which defaults to Euclidian distance) are also
# close in ordering. More details can be found in the :ref:`User Guide<OPTICS>`
#
# These ideas are better visualized by looking at a reachability plot. This is
# obtained by plotting the ordering of the samples on the horizontal axis
# against the corresponding reachability distance of the sample on the vertical
# axis. We will also plot horizontal lines of the form `y=eps` for the values
# of `eps` used previously. The reason for this will be made clear very soon.

_, ax = plt.subplots(figsize=(10, 4))

space = np.arange(len(X))
reachability = optics.reachability_[optics.ordering_]

labels = optics.labels_[optics.ordering_]

colors = cluster_colours(labels)

for k, col in zip(set(labels), colors):
    if k == -1:
        continue

    ax.plot(space[labels == k], reachability[labels == k], color=tuple(col), alpha=0.7)

for eps in eps_values:
    ax.plot(space, np.full_like(space, eps, dtype=float), "k-.", alpha=0.5)

ax.set_title("Reachability Plot")
plt.tight_layout()

# %%
# Observe that points in the two dense clusters have low reachability
# distances, because points in a dense cluster are close to core
# points in that cluster. On the other hand, samples in the sparse
# clusters have higher reachability distances.
#
# Clusters can be extracted from this plot in two ways, either Xi or DBSCAN
# clustering. The user can choose which method by setting the `cluster_method`
# parameter to either `xi` (default) or `dbscan`.
#
# The Xi clustering method uses the steep slopes within the reachability
# plot to determine which points are in the same cluster. One can specify
# this by setting the `xi` parameter. More details on this parameter
# may be found the :ref:`User Guide <OPTICS>`.
#
# Extracting DBSCAN clusters
# ++++++++++++++++++++++++++
# We can also extract clustering from the reachability plot to produce results
# similar to `DBSCAN` for a given value of `eps`. This can be done manually by
# plotting the horizontal line `y=eps` on the reachability plot and reading the
# plot left to right. Breaks in the graph signify new clusters, and points
# above the line `y=eps` are treated as noise.
#
# We may also do this computationally using `cluster_optics_dbscan()` or
# setting `cluster_method='dbscan'` and an appropriate value for `eps` and
# other hyperparameters before model fitting. We will discuss only
# `cluster_optics_dbscan()`. This function takes the reachability, ordering,
# and core distances produced from running OPTICS, and an `eps` parameter
# analogous to its counterpart in DBSCAN. It produces a clustering with results
# similar if one were to run DBSCAN on the same dataset with similar settings
# for `eps` and `min_samples`. The runtime is efficient, it is linear in the
# number of samples whereas running DBSCAN from scratch is quadratic in its
# worst case.
#
# We run `cluster_optics_dbscan` over several `eps` values. For the following
# plots, each corresponding to a specific value of `eps`, observe the
# similarity to the corresponding clustering we previously we made using
# DBSCAN. Also consider how the clustering can be obtained from the
# reachability plot with the threshold `y=eps`.

fig, axes = plt.subplots(len(eps_values), 1, figsize=(10, 12))


def plot_optics_dbscan(optic_clustering, eps, ax):
    label_eps = cluster_optics_dbscan(
        reachability=optics.reachability_,
        core_distances=optics.core_distances_,
        ordering=optics.ordering_,
        eps=eps,
    )
    plot(X, label_eps, parameters={"eps": eps}, ax=ax)


for idx, eps in enumerate(eps_values):
    plot_optics_dbscan(optics, eps, axes[idx])
