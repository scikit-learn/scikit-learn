"""
=======================================
Clustering text documents using k-means
=======================================

This is an example showing how the scikit-learn API can be used to cluster
documents by topics using a `Bag of Words approach
<https://en.wikipedia.org/wiki/Bag-of-words_model>`_.

Two algorithms are demoed: :func:`~sklearn.cluster.KMeans` and its
more scalable variant, :func:`~sklearn.cluster.MiniBatchKMeans`. Additionally,
latent semantic analysis is used to reduce dimensionality and discover latent
patterns in the data.

This example uses two different text vectorizers: a
:func:`~sklearn.feature_extraction.text.TfidfVectorizer` and a
:func:`~sklearn.feature_extraction.text.HashingVectorizer`. See the example
notebook :ref:`sphx_glr_auto_examples_text_plot_hashing_vs_dict_vectorizer.py`
for more information on vectorizers and a comparison of their processing times.

For document analysis via a supervised learning approach, see the example script
:ref:`sphx_glr_auto_examples_text_plot_document_classification_20newsgroups.py`.

"""

# Author: Peter Prettenhofer <peter.prettenhofer@gmail.com>
#         Lars Buitinck
#         Arturo Amor <david-arturo.amor-quiroz@inria.fr>
# License: BSD 3 clause

# %%
# Load Data
# =========
#
# We load data from :ref:`20newsgroups_dataset`, which comprises around 18,000
# newsgroups posts on 20 topics. For illustrative purposes and to reduce the
# computational cost, we select a subset of 4 topics only accounting for around
# 3,400 documents. See the example
# :ref:`sphx_glr_auto_examples_text_plot_document_classification_20newsgroups.py`
# to gain intuition on the overlap of such topics.
#
# Notice that, by default, the text samples contain some message metadata such
# as `"headers"`, `"footers"` (signatures) and `"quotes"` to other posts. The
# :func:`~sklearn.datasets.fetch_20newsgroups` function therefore accepts a
# parameter named `remove` to attempt stripping such information that can make
# the clustering problem "too easy".

import numpy as np
from sklearn.datasets import fetch_20newsgroups

categories = [
    "alt.atheism",
    "talk.religion.misc",
    "comp.graphics",
    "sci.space",
]

dataset = fetch_20newsgroups(
    remove=("headers", "footers", "quotes"),
    subset="all",
    categories=categories,
    shuffle=True,
    random_state=42,
)

labels = dataset.target
true_k = np.unique(labels).shape[0]

print(f"{len(dataset.data)} documents - {true_k} categories")

# %%
# Comparing Clustering
# ====================
#
# Here we define a function to score different clustering pipelines
# using several metrics.
#
# Given the knowledge of the ground truth class assignments of the samples, it
# is possible to define metrics to test for different properties of a cluster:
#
# - homogeneity, which quantifies cluster containing only members of a single
# class;
#
# - completeness, which quantifies cluster having members of a given class are
# assigned to the same cluster;
#
# - V-measure, the harmonic mean of completeness and homogeneity;
#
# - Rand-Index, which measures the similarity between pairs of clusters
#
# - Adjusted Rand-Index, a change-adjusted Rand-Index
#
# If the ground truth labels are not known, evaluation can only be performed
# using the model results itself. In that case, the Silhouette Coefficient
# comes in handy.
#
# For more reference, see :ref:`_clustering_evaluation`.

from collections import defaultdict
from sklearn import metrics
from time import time

scores = defaultdict(list)
train_times = defaultdict(list)


def fit_and_evaluate(km, X, name=None):
    t0 = time()
    km.fit(X)
    train_time = time() - t0
    name = km.__class__.__name__ if name is None else name
    train_times["estimator"].append(name)
    train_times["train time"].append(train_time)
    scores["estimator"].append(name)
    scores["Homogeneity"].append(metrics.homogeneity_score(labels, km.labels_))
    scores["Completeness"].append(metrics.completeness_score(labels, km.labels_))
    scores["V-measure"].append(metrics.v_measure_score(labels, km.labels_))
    scores["Adjusted Rand-Index"].append(
        metrics.adjusted_rand_score(labels, km.labels_)
    )
    scores["Silhouette Coefficient"].append(
        metrics.silhouette_score(X, km.labels_, sample_size=2000)
    )
    return


# %%
# Feature Extraction using sparse vectorizers
# ===========================================
#
# Two feature extraction methods are used in this example:
#
# - :func:`~sklearn.feature_extraction.text.TfidfVectorizer` uses a in-memory
#   vocabulary (a python dict) to map the most frequent words to features
#   indices and hence compute a word occurrence frequency (sparse) matrix. The
#   word frequencies are then reweighted using the Inverse Document Frequency
#   (IDF) vector collected feature-wise over the corpus.
#
# - :func:`~sklearn.feature_extraction.text.HashingVectorizer` hashes word
#   occurrences to a fixed dimensional space, possibly with collisions. The word
#   count vectors are then normalized to each have l2-norm equal to one
#   (projected to the euclidean unit-sphere) which seems to be important for
#   k-means to work in high dimensional space.
#
# TfidfVectorizer
# ---------------
#
# We first benchmark the estimators using a dictionary vectorizer along with an
# IDF normalization as provided by
# :func:`~sklearn.feature_extraction.text.TfidfVectorizer`

from sklearn.feature_extraction.text import TfidfVectorizer

vectorizer = TfidfVectorizer(
    max_df=0.5,
    min_df=2,
    stop_words="english",
)
t0 = time()
X = vectorizer.fit_transform(dataset.data)

print(f"vectorizing done in {time() - t0:.3f} s")
print(f"n_samples: {X.shape[0]}, n_features: {X.shape[1]}")

# %%
# **Clustering sparse data with MiniBatchKMeans**
#
# As both :func:`~sklearn.cluster.KMeans` and
# :func:`~sklearn.cluster.MiniBatchKMeans` optimize a non-convex objective
# function, their clustering will likely be the a sub-optimal one.
# Several runs with independent random initiations are performed using the
# `n_init` and the clustering with the smallest inertia is chosen.

from sklearn.cluster import MiniBatchKMeans

km = MiniBatchKMeans(
    n_clusters=true_k,
    init="k-means++",
    n_init=5,
    init_size=1000,
    batch_size=1000,
    random_state=0,
)

fit_and_evaluate(km, X)

# %%
# **Clustering sparse data with KMeans**

from sklearn.cluster import KMeans

km = KMeans(
    n_clusters=true_k,
    init="k-means++",
    max_iter=100,
    n_init=5,
    random_state=0,
)

fit_and_evaluate(km, X)

# %%
# **Performing dimensionality reduction using LSA**
#
# Since LSA/SVD results are not normalized, we redo the normalization to improve
# the :func:`~sklearn.cluster.KMeans` result.

from sklearn.decomposition import TruncatedSVD
from sklearn.pipeline import make_pipeline
from sklearn.preprocessing import Normalizer

svd = TruncatedSVD(n_components=100)
lsa = make_pipeline(svd, Normalizer(copy=False))
X = lsa.fit_transform(X)
fit_and_evaluate(km, X, name=km.__class__.__name__ + " with LSA")

explained_variance = svd.explained_variance_ratio_.sum()
print(f"Explained variance of the SVD step: {explained_variance * 100:.1f}%")

# %%
# **Top terms per cluster**

original_space_centroids = svd.inverse_transform(km.cluster_centers_)
order_centroids = original_space_centroids.argsort()[:, ::-1]
terms = vectorizer.get_feature_names_out()

for i in range(true_k):
    print("Cluster %d:" % i, end="")
    for ind in order_centroids[i, :10]:
        print(" %s" % terms[ind], end="")
    print()

# %%
# HashingVectorizer
# -----------------
# A similar experiment can be done using a
# :func:`~sklearn.feature_extraction.text.HashingVectorizer` instance, which
# does not provide IDF weighting as this is a stateless model (the fit method
# does nothing). Here we illustrate the effect of Tfidf weighting by
# benchmarking the function :func:`~sklearn.cluster.KMeans` with and without
# such normalization.

from sklearn.feature_extraction.text import HashingVectorizer

hasher = HashingVectorizer(
    stop_words="english",
    alternate_sign=False,
    norm=None,
)

t0 = time()
X = hasher.fit_transform(dataset.data)
print(f"vectorizing done in {time() - t0:.3f} s")
print(f"n_samples: {X.shape[0]}, n_features: {X.shape[1]}")

fit_and_evaluate(km, X, name=km.__class__.__name__ + " with\nsimple hashing")

# %%
# When IDF weighting is needed it can be added by pipelining the
# :func:`~sklearn.feature_extraction.text.HashingVectorizer` output to a
# :func:`~sklearn.feature_extraction.text.TfidfTransformer` instance.

from sklearn.feature_extraction.text import TfidfTransformer

vectorizer = make_pipeline(hasher, TfidfTransformer())

t0 = time()
X = vectorizer.fit_transform(dataset.data)
print(f"vectorizing done in {time() - t0:.3f} s")

fit_and_evaluate(km, X, name=km.__class__.__name__ + " with\nTfidf-scaled hashing")

# %%
# Plot unsupervised evaluation metrics
# ====================================

import matplotlib
import matplotlib.pyplot as plt
import pandas as pd
from itertools import cycle

df = pd.DataFrame(scores).set_index("estimator")

ax = df.plot.barh()
ax.set_ylabel("")
_ = ax.set_xlabel("Clustering scores")

# %%
# It can be noticed that :func:`~sklearn.cluster.KMeans` (and
# :func:`~sklearn.cluster.MiniBatchKMeans`) are very sensitive to feature
# scaling and that in this case the IDF weighting helps improve the quality of
# the clustering by quite a lot as measured against the "ground truth" provided
# by the class label assignments of :ref:`20newsgroups_dataset`.
#
# This improvement is not visible in the Silhouette Coefficient which is small
# for both as this measure seem to suffer from the phenomenon called
# "Concentration of Measure" or "Curse of Dimensionality" for high dimensional
# datasets such as text data. Other measures such as V-measure and Adjusted Rand
# Index are information-theory-based evaluation scores: as they are only based
# on cluster assignments rather than distances, hence not affected by the curse
# of dimensionality.


# %%
# Plot processing time
# ====================

fig, ax = plt.subplots(figsize=(12, 6))

y_pos = np.arange(len(train_times["estimator"]))
ax.barh(y_pos, train_times["train time"], align="center")
ax.set_yticks(y_pos)
ax.set_yticklabels(train_times["estimator"])
ax.invert_yaxis()
_ = ax.set_xlabel("processing times (s)")

# %%
# :func:`~sklearn.cluster.MiniBatchKMeans` requires less time since
# it runs computations on batches of data. This comes at the expense
# of clustering quality.
