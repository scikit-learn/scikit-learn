"""
=======================================
Clustering text documents using k-means
=======================================

This is an example showing how the scikit-learn API can be used to cluster
documents by topics using a `Bag of Words approach
<https://en.wikipedia.org/wiki/Bag-of-words_model>`_.

Two algorithms are demoed: ordinary :func:`~sklearn.cluster.KMeans` and its more
scalable cousin :func:`~sklearn.cluster.MiniBatchKMeans`. Additionally, latent
semantic analysis is used to reduce dimensionality and discover latent patterns
in the data.

This example uses two different text vectorizers: a
:func:`~sklearn.feature_extraction.text.TfidfVectorizer` and a
:func:`~sklearn.feature_extraction.text.HashingVectorizer`. See the example
notebook :ref:`sphx_glr_auto_examples_text_plot_hashing_vs_dict_vectorizer.py`.
For more information on vectorizers and a comparison of their processing times.

Clustering is unsupervised learning. For the supervised learning version of
document analysis, see the example script
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
# We load data from :ref:`20newsgroups_dataset`, which comprises around 18000
# newsgroups posts on 20 topics. For illustrative purposes and reducing the
# computational cost, we select a subset of 4 topics only. See the example
# :ref:`sphx_glr_auto_examples_text_plot_document_classification_20newsgroups.py`
# to gain intuition on the overlap of such topics.
#
# Notice that, by default, the text samples contain some message metadata such
# as `"headers"`, `"footers"` (signatures) and `"quotes"` to other posts. The
# fetch_20newsgroups function therefore accepts a parameter named `remove` to
# attempt stripping such information that can make the clustering problem “too
# easy”.

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
# Benchmark Clustering
# ====================
#
# Here we define a function to score different clustering pipelines. For such
# purpose we evaluate several metrics.
#
# Given the knowledge of the ground truth class assignments of the samples, it
# is possible to define metrics to test for different properties of a cluster:
#
# - homogeneity: each cluster contains only members of a single class;
#
# - completeness: all members of a given class are assigned to the same cluster;
#
# - V-measure: the harmonic mean of completeness and homogeneity;
#
# - Rand-Index: similarity between pairs of clusters
#   `RI = (number of agreeing pairs) / (number of pairs)` ;
#
# - Adjusted Rand-Index: Rand-Index adjusted for chance
#   `ARI = (RI - Expected_RI) / (max(RI) - Expected_RI)` .
#
# If the ground truth labels are not known, evaluation must be performed using
# the model itself. The Silhouette Coefficient is an example of such metrics,
# where the score is higher when clusters are dense and well separated.

from collections import defaultdict
from sklearn import metrics

scores = defaultdict(list)
train_times = defaultdict(list)


def benchmark(km, X, custom_name=False):
    t0 = time()
    km.fit(X)
    train_time = time() - t0
    if custom_name:
        scores["estimator"].append(km.__class__.__name__ + custom_name)
        train_times["estimator"].append(km.__class__.__name__ + custom_name)
    else:
        scores["estimator"].append(km.__class__.__name__)
        train_times["estimator"].append(km.__class__.__name__)
    train_times["train time"].append(train_time)
    scores["Homogeneity"].append(metrics.homogeneity_score(labels, km.labels_))
    scores["Completeness"].append(metrics.completeness_score(labels, km.labels_))
    scores["V-measure"].append(metrics.v_measure_score(labels, km.labels_))
    scores["Adjusted Rand-Index"].append(
        metrics.adjusted_rand_score(labels, km.labels_)
    )
    scores["Silhouette Coefficient"].append(
        metrics.silhouette_score(X, km.labels_, sample_size=1000)
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
#   (projected to the euclidean unit-ball) which seems to be important for
#   k-means to work in high dimensional space.
#
#   :func:`~sklearn.feature_extraction.text.HashingVectorizer` does not provide
#   IDF weighting as this is a stateless model (the fit method does nothing).
#   When IDF weighting is needed it can be added by pipelining its output to a
#   :func:`~sklearn.feature_extraction.text.TfidfTransformer` instance.
#
# TfidfVectorizer
# ---------------
#
# We first benchmark the estimators using a dictionary vectorizer along with an
# IDF normalization as provided by
# :func:`~sklearn.feature_extraction.text.TfidfVectorizer`

from sklearn.feature_extraction.text import TfidfVectorizer
from time import time

t0 = time()

vectorizer = TfidfVectorizer(
    max_df=0.5,
    min_df=2,
    stop_words="english",
)
X = vectorizer.fit_transform(dataset.data)

print(f"done in {time() - t0}s")
print(f"n_samples: {X.shape[0]}, n_features: {X.shape[1]}")

# %%
# **Clustering sparse data with MiniBatchKMeans**
#
# As both KMeans and MiniBatchKMeans optimize a non-convex objective function,
# they will likely end up in a local optimum. Several runs with independent
# random init might be necessary to get a good convergence. The parameter
# `n_init` controls the number of such initializations.

from sklearn.cluster import MiniBatchKMeans

km = MiniBatchKMeans(
    n_clusters=true_k,
    init="k-means++",
    n_init=5,
    init_size=1000,
    batch_size=1000,
)

benchmark(km, X)

# %%
# **Clustering sparse data with KMeans**

from sklearn.cluster import KMeans

km = KMeans(
    n_clusters=true_k,
    init="k-means++",
    max_iter=100,
    n_init=5,
)

benchmark(km, X)

# %%
# **Performing dimensionality reduction using LSA**
#
# Since LSA/SVD results are not normalized, we redo the normalization to improve
# the k-means result.

from sklearn.decomposition import TruncatedSVD
from sklearn.pipeline import make_pipeline
from sklearn.preprocessing import Normalizer

svd = TruncatedSVD(n_components=100)
normalizer = Normalizer(copy=False)
lsa = make_pipeline(svd, normalizer)
X = lsa.fit_transform(X)
benchmark(km, X, custom_name=" with LSA")

explained_variance = svd.explained_variance_ratio_.sum()
print(f"Explained variance of the SVD step: {int(explained_variance * 100)}%")

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
# :func:`~sklearn.feature_extraction.text.HashingVectorizer` instance. Here
# we illustrate the effect of Tfidf weighting by benchmarking the KMeans
# algorithm with and without such normalization.

from sklearn.feature_extraction.text import HashingVectorizer
from sklearn.feature_extraction.text import TfidfTransformer

hasher = HashingVectorizer(
    stop_words="english",
    alternate_sign=False,
    norm=None,
)
vectorizer = make_pipeline(hasher, TfidfTransformer())

X = hasher.fit_transform(dataset.data)
benchmark(km, X, custom_name=" with\nsimple hashing")

X = vectorizer.fit_transform(dataset.data)
benchmark(km, X, custom_name=" with\nhashed Tfidf")

# %%
# Plot unsupervised evaluation metrics
# ====================================

import matplotlib
import matplotlib.pyplot as plt
import pandas as pd
from itertools import cycle

df = pd.DataFrame(scores)
df = df.set_index("estimator")

bar_size = 0.25
padding = 0.75
bars_per_group = len(df.columns)
cmap = matplotlib.cm.get_cmap("rainbow")
y_locs = np.arange(len(df.index)) * (bars_per_group * bar_size + padding)
y_array = [y_locs + (j - bars_per_group / 2) * bar_size for j in range(bars_per_group)]
colours = cycle([cmap(index) for index in np.linspace(0, 1, bars_per_group)])

fig, ax = plt.subplots(figsize=(10, 8))
for i, estimator in enumerate(df.index):
    for j, label in enumerate(df.columns):
        ax.barh(
            y_array[j][i],
            df.loc[estimator][label],
            height=bar_size,
            color=next(colours),
        )
ax.set(
    yticks=y_locs - bar_size / 2,
    yticklabels=df.index,
    ylim=[
        0 - bars_per_group * bar_size,
        len(df.index) * (bars_per_group * bar_size + padding)
        - bars_per_group * bar_size,
    ],
)
ax.legend(loc="lower right", labels=df.columns)
_ = ax.set_title("Clustering scores")

# %%
# It can be noted that k-means (and minibatch k-means) are very sensitive to
# feature scaling and that in this case the IDF weighting helps improve the
# quality of the clustering by quite a lot as measured against the "ground
# truth" provided by the class label assignments of
# :ref:`20newsgroups_dataset`.
#
# This improvement is not visible in the Silhouette Coefficient which is small
# for both as this measure seem to suffer from the phenomenon called
# "Concentration of Measure" or "Curse of Dimensionality" for high dimensional
# datasets such as text data. Other measures such as V-measure and Adjusted Rand
# Index are information-theory-based evaluation scores: as they are only
# based on cluster assignments rather than distances, hence not affected by the
# curse of dimensionality.


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
# :func:`~sklearn.cluster.MiniBatchKMeans` requires less time due to scalability.
# That is, at the expense of clustering quality.
