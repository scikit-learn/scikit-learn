"""
==========================================================
Adjustment for chance in clustering performance evaluation
==========================================================

This notebook explores the impact of random labeling on the behavior of some
clustering evaluation metrics when we vary the number of clusters in:
 - the labeling to evaluate;
 - the reference labeling we compare to.

"""

# Author: Olivier Grisel <olivier.grisel@ensta.org>
# License: BSD 3 clause

# %%
# Quantifying the score of random labeling
# ----------------------------------------
#
# In this section we define a function that uses several metrics to score 2
# uniformly-distributed random labelings.
#
# Both random labelings have the same number of clusters for each value possible
# value in `n_clusters_range`. When fixed_n_classes is not `None`, the first
# labeling is considered a ground truth class assignment with fixed number of
# classes.

import numpy as np


def uniform_labelings_scores(
    score_func, n_samples, n_clusters_range, n_runs=5, seed=42
):
    random_labels = np.random.RandomState(seed).randint
    scores = np.zeros((len(n_clusters_range), n_runs))

    for i, k in enumerate(n_clusters_range):
        for j in range(n_runs):
            labels_a = random_labels(low=0, high=k, size=n_samples)
            labels_b = random_labels(low=0, high=k, size=n_samples)
            scores[i, j] = score_func(labels_a, labels_b)
    return scores


# %%
# Clustering algorithms are fundamentally unsupervised learning methods.
# However, since we assigned class labels for this synthetic dataset, it is
# possible to use evaluation metrics that leverage this "supervised" ground
# truth information to quantify the quality of the resulting clusters. Examples
# of such metrics are the following:
#
# - V-measure, the harmonic mean of completeness and homogeneity;
#
# - Rand-Index, which measures how frequently pairs of data points are grouped
#   consistently according to the result of the clustering algorithm and the
#   ground truth class assignment;
#
# - Adjusted Rand-Index (ARI), a chance-adjusted Rand-Index such that random
#   cluster assignment has an ARI of 0.0 in expectation;
#
# - Mutual Information (MI), determines how different the joint distribution of
#   the pair `(X,Y)` is from the product of the marginal distributions of `X`
#   and `Y`;
#
# - Adjusted Mutual Information (AMI), a chance-adjusted Mutual Information.
#   Similarly to ARI, random cluster assignment has an AMI of 0.0 in
#   expectation;
#
# - Normalized Mutual Information (NMI), a Mutual Information defined between 0
#   (no mutual information) and 1 (perfect correlation). It is not adjusted for
#   chance.
#
# For more information, see the :ref:`clustering_evaluation` module.

from sklearn import metrics

score_funcs = [
    ("V-measure", metrics.v_measure_score),
    ("Rand-Index", metrics.rand_score),
    ("ARI", metrics.adjusted_rand_score),
    ("MI", metrics.mutual_info_score),
    ("AMI", metrics.adjusted_mutual_info_score),
    ("NMI", metrics.normalized_mutual_info_score),
]

# %%
# Plot clustering scores
# ----------------------
#
# We first score 2 independent random clusterings with equal cluster number

import matplotlib.pyplot as plt
from time import time

n_samples = 100
n_clusters_range = np.linspace(2, n_samples, 10).astype(int)

plt.figure(1)

plots = []
names = []
for score_name, score_func in score_funcs:
    print(
        "Computing %s for %d values of n_clusters and n_samples=%d"
        % (score_name, len(n_clusters_range), n_samples)
    )

    t0 = time()
    scores = uniform_labelings_scores(score_func, n_samples, n_clusters_range)
    print("done in %0.3fs" % (time() - t0))
    plots.append(
        plt.errorbar(
            n_clusters_range, np.median(scores, axis=1), scores.std(axis=1), alpha=0.8
        )[0]
    )
    names.append(score_name)

plt.title(
    "Clustering measures for 2 random uniform labelings\nwith equal number of clusters"
)
plt.xlabel("Number of clusters (Number of samples is fixed to %d)" % n_samples)
plt.ylabel("Score value")
plt.legend(plots, names)
plt.ylim(bottom=-0.05, top=1.05)
plt.show()

# %%
# Non-adjusted measures such as the V-Measure show a dependency between the
# number of clusters and the number of samples: the mean V-Measure of random
# labeling increases significantly as the number of clusters is closer to the
# total number of samples used to compute the measure.
#
# Adjusted for chance measure, such as ARI, display some random variations
# centered around a mean score of 0.0 for any number of samples and clusters.
#
# We can also score a random labeling with varying n_clusters against the ground
# class labels with a fixed number of samples (1000).

n_samples = 1000
n_clusters_range = np.linspace(2, 100, 10).astype(int)
n_classes = 10


def fixed_classes_uniform_labelings_scores(
    score_func, n_samples, n_clusters_range, n_classes=10, n_runs=5, seed=42
):
    random_labels = np.random.RandomState(seed).randint
    scores = np.zeros((len(n_clusters_range), n_runs))
    labels_a = random_labels(low=0, high=n_classes, size=n_samples)

    for i, k in enumerate(n_clusters_range):
        for j in range(n_runs):
            labels_b = random_labels(low=0, high=k, size=n_samples)
            scores[i, j] = score_func(labels_a, labels_b)
    return scores


plt.figure(2)

plots = []
names = []
for score_name, score_func in score_funcs:
    print(
        "Computing %s for %d values of n_clusters and n_samples=%d"
        % (score_name, len(n_clusters_range), n_samples)
    )

    t0 = time()
    scores = fixed_classes_uniform_labelings_scores(
        score_func, n_samples, n_clusters_range, n_classes=n_classes
    )
    print("done in %0.3fs" % (time() - t0))
    plots.append(
        plt.errorbar(
            n_clusters_range, scores.mean(axis=1), scores.std(axis=1), alpha=0.8
        )[0]
    )
    names.append(score_name)

plt.title(
    "Clustering measures for random uniform labeling\n"
    "against reference assignment with %d classes" % n_classes
)
plt.xlabel("Number of clusters (Number of samples is fixed to %d)" % n_samples)
plt.ylabel("Score value")
plt.ylim(bottom=-0.05, top=1.05)
plt.legend(plots, names)
plt.show()

# %%
# We observe similar results as for the first experiment: adjusted for chance
# metrics stay constantly near zero while other metrics tend to get larger with
# finer-grained labelings.
#
# Only adjusted measures can hence safely be used as a consensus index to
# evaluate the average stability of clustering algorithms for a given value of k
# on various overlapping sub-samples of the dataset.
#
# Non-adjusted clustering evaluation metric can therefore be misleading as they
# output large values for fine-grained labelings, one could be lead to think
# that the labeling has captured meaningful groups while they can be totally
# random. In particular, such non-adjusted metrics should not be used to compare
# the results of different clustering algorithms that output a different number
# of clusters.
