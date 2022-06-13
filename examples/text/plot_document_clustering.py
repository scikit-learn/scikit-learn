"""
=======================================
Clustering text documents using k-means
=======================================

This is an example showing how the scikit-learn API can be used to cluster
documents by topics using a `Bag of Words approach
<https://en.wikipedia.org/wiki/Bag-of-words_model>`_.

Two algorithms are demoed: :class:`~sklearn.cluster.KMeans` and its
more scalable variant, :class:`~sklearn.cluster.MiniBatchKMeans`. Additionally,
latent semantic analysis is used to reduce dimensionality and discover latent
patterns in the data.

This example uses two different text vectorizers: a
:class:`~sklearn.feature_extraction.text.TfidfVectorizer` and a
:class:`~sklearn.feature_extraction.text.HashingVectorizer`. See the example
notebook :ref:`sphx_glr_auto_examples_text_plot_hashing_vs_dict_vectorizer.py`
for more information on vectorizers and a comparison of their processing times.

For document analysis via a supervised learning approach, see the example script
:ref:`sphx_glr_auto_examples_text_plot_document_classification_20newsgroups.py`.

"""

# Author: Peter Prettenhofer <peter.prettenhofer@gmail.com>
#         Lars Buitinck
#         Olivier Grisel <olivier.grisel@ensta.org>
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
# as `"headers"`, `"footers"` (signatures) and `"quotes"` to other posts. We use
# the `remove` parameter from :func:`~sklearn.datasets.fetch_20newsgroups` to
# strip those features and have a more sensible clustering problem.

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
# In this section we define a function to score different clustering pipelines
# using several metrics.
#
# Clustering algorithms are fundamentally unsupervised learning methods.
# However, since we happen to have class labels for this specific dataset, it is
# possible to use evaluation metrics that leverage this "supervised" ground
# truth information to quantify the quality of the resulting clusters. Examples
# of such metrics are the following:
#
# - homogeneity, which quantifies cluster containing only members of a single
#   class;
#
# - completeness, which quantifies cluster having members of a given class are
#   assigned to the same cluster;
#
# - V-measure, the harmonic mean of completeness and homogeneity;
#
# - Rand-Index, which measures the similarity between pairs of clusters;
#
# - Adjusted Rand-Index, a change-adjusted Rand-Index.
#
# If the ground truth labels are not known, evaluation can only be performed
# using the model results itself. In that case, the Silhouette Coefficient
# comes in handy.
#
# For more reference, see :ref:`clustering_evaluation`.

from collections import defaultdict
from sklearn import metrics
from time import time


evaluations = []
evaluations_std = []


def fit_and_evaluate(km, X, name=None, n_runs=5):
    name = km.__class__.__name__ if name is None else name

    train_times = []
    scores = defaultdict(list)
    for seed in range(n_runs):
        km.set_params(random_state=seed)
        t0 = time()
        km.fit(X)
        train_times.append(time() - t0)
        scores["Homogeneity"].append(metrics.homogeneity_score(labels, km.labels_))
        scores["Completeness"].append(metrics.completeness_score(labels, km.labels_))
        scores["V-measure"].append(metrics.v_measure_score(labels, km.labels_))
        scores["Adjusted Rand-Index"].append(
            metrics.adjusted_rand_score(labels, km.labels_)
        )
        scores["Silhouette Coefficient"].append(
            metrics.silhouette_score(X, km.labels_, sample_size=2000)
        )
    train_times = np.asarray(train_times)

    print(f"clustering done in {train_times.mean():.1f} ± {train_times.std():.1f} s ")
    evaluation = {
        "estimator": name,
        "train_time": train_times.mean(),
    }
    evaluation_std = {
        "estimator": name,
        "train_time": train_times.std(),
    }
    for score_name, score_values in scores.items():
        mean_score, std_score = np.mean(score_values), np.std(score_values)
        print(f"{score_name}: {mean_score:.3f} ± {std_score:.3f}")
        evaluation[score_name] = mean_score
        evaluation_std[score_name] = std_score
    evaluations.append(evaluation)
    evaluations_std.append(evaluation_std)


# %%
# Feature Extraction using sparse vectorizers
# ===========================================
#
# Two feature extraction methods are used in this example:
#
# - :class:`~sklearn.feature_extraction.text.TfidfVectorizer` uses an in-memory
#   vocabulary (a python dict) to map the most frequent words to features
#   indices and hence compute a word occurrence frequency (sparse) matrix. The
#   word frequencies are then reweighted using the Inverse Document Frequency
#   (IDF) vector collected feature-wise over the corpus.
#
# - :class:`~sklearn.feature_extraction.text.HashingVectorizer` hashes word
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
# :class:`~sklearn.feature_extraction.text.TfidfVectorizer`.

from sklearn.feature_extraction.text import TfidfVectorizer

vectorizer = TfidfVectorizer(
    max_df=0.5,
    min_df=5,
    stop_words="english",
)
t0 = time()
X_tfidf = vectorizer.fit_transform(dataset.data)

print(f"vectorizing done in {time() - t0:.3f} s")
print(f"n_samples: {X_tfidf.shape[0]}, n_features: {X_tfidf.shape[1]}")

# %%
# **Clustering sparse data with KMeans**
#
# As both :class:`~sklearn.cluster.KMeans` and
# :class:`~sklearn.cluster.MiniBatchKMeans` optimize a non-convex objective
# function, their clustering will likely be sub-optimal.
# Several runs with independent random initiations are performed using the
# `n_init` and the clustering with the smallest inertia is chosen.
#
# Here we use a `"random"` init instead of the smart `"k-means++"` init. The
# `"k-means++"` init picks centroids with proba inverse to square distance,
# which means it will pick points extremely isolated if there are any, which
# then stay their own centroids all along. That is usually the case for text
# clustering, where the vectorized space is high dimensional and very sparse.

from sklearn.cluster import KMeans

kmeans = KMeans(
    n_clusters=true_k,
    init="random",
    max_iter=100,
    n_init=5,
)

fit_and_evaluate(kmeans, X_tfidf, name="KMeans(init='random')\non tf-idf vectors")

# %%
# **Performing dimensionality reduction using LSA**
#
# The `"k-means++"` init can still be used as long as the dimension of the
# vectorized space is reduced first. For such purpose we use
# :class:`~sklearn.decomposition.TruncatedSVD`, which works on term count/tf-idf
# matrices. Since SVD results are not normalized, we redo the normalization to
# improve the :class:`~sklearn.cluster.KMeans` result. Such process is known as
# latent semantic analysis (LSA).

from sklearn.decomposition import TruncatedSVD
from sklearn.pipeline import make_pipeline
from sklearn.preprocessing import Normalizer


lsa = make_pipeline(TruncatedSVD(n_components=100), Normalizer(copy=False))
t0 = time()
X_lsa = lsa.fit_transform(X_tfidf)
explained_variance = lsa[0].explained_variance_ratio_.sum()

print(f"LSA done in {time() - t0:.3f} s")
print(f"Explained variance of the SVD step: {explained_variance * 100:.1f}%")

# %%
# We are now ready to use the `"k-means++"` init. As it requires a single
# initialization, the processing time will be reduced for both
# :class:`~sklearn.cluster.KMeans` and :class:`~sklearn.cluster.MiniBatchKMeans`.

kmeans = KMeans(
    n_clusters=true_k,
    init="k-means++",
    max_iter=100,
    n_init=1,
)

fit_and_evaluate(
    kmeans, X_lsa, name="KMeans(init='kmeans++')\nwith LSA on tf-idf vectors"
)

# %%
# We repeat the experiment with :class:`~sklearn.cluster.MiniBatchKMeans`.

from sklearn.cluster import MiniBatchKMeans

minibatch_kmeans = MiniBatchKMeans(
    n_clusters=true_k,
    init="k-means++",
    n_init=1,
    init_size=1000,
    batch_size=1000,
)

fit_and_evaluate(
    minibatch_kmeans,
    X_lsa,
    name="MiniBatchKMeans(init='kmeans++')\nwith LSA on tf-idf vectors",
)

# %%
# **Top terms per cluster**
#
# Since :class:`~sklearn.feature_extraction.text.TfidfVectorizer` and can be
# inverted we can identify the cluster centers, which provide an intuition of
# the most predictive words **per cluster**. See the example script
# :ref:`sphx_glr_auto_examples_text_plot_document_classification_20newsgroups.py`
# for a comparison with the most predictive words **per class**.

original_space_centroids = lsa[0].inverse_transform(kmeans.cluster_centers_)
order_centroids = original_space_centroids.argsort()[:, ::-1]
terms = vectorizer.get_feature_names_out()

for i in range(true_k):
    print(f"Cluster {i}: ", end="")
    for ind in order_centroids[i, :10]:
        print(f"{terms[ind]} ", end="")
    print()

# %%
# HashingVectorizer
# -----------------
# A similar experiment can be done using a
# :class:`~sklearn.feature_extraction.text.HashingVectorizer` instance, which
# does not provide IDF weighting as this is a stateless model (the fit method
# does nothing). When IDF weighting is needed it can be added by pipelining the
# :class:`~sklearn.feature_extraction.text.HashingVectorizer` output to a
# :class:`~sklearn.feature_extraction.text.TfidfTransformer` instance. In this
# case we also add LSA to the pipeline to evaluate the `"k-means++"` init.

from sklearn.feature_extraction.text import HashingVectorizer
from sklearn.feature_extraction.text import TfidfTransformer

hasher = HashingVectorizer(
    stop_words="english",
    n_features=50_000,
)

lsa_vectorizer = make_pipeline(
    hasher,
    TfidfTransformer(),
    TruncatedSVD(n_components=100, random_state=0),
    Normalizer(copy=False),
)

t0 = time()
X_hashed_lsa = lsa_vectorizer.fit_transform(dataset.data)
print(f"vectorizing done in {time() - t0:.3f} s")

fit_and_evaluate(
    kmeans, X_hashed_lsa, name="KMeans(init='kmeans++')\nwith LSA on hashed vectors"
)

# %%
fit_and_evaluate(
    minibatch_kmeans,
    X_hashed_lsa,
    name="MiniBatchKMeans(init='kmeans++')\nwith LSA on hashed vectors",
)

# %%
# Clustering evaluation summary
# ==============================

import pandas as pd
import matplotlib.pyplot as plt

fig, (ax0, ax1) = plt.subplots(ncols=2, figsize=(16, 6), sharey=True)

df = pd.DataFrame(evaluations).set_index("estimator")
df_std = pd.DataFrame(evaluations_std).set_index("estimator")

df.drop(
    ["train_time"],
    axis="columns",
).plot.barh(ax=ax0, xerr=df_std)
ax0.set_xlabel("Clustering scores")
ax0.set_ylabel("")

df["train_time"].plot.barh(ax=ax1, xerr=df_std["train_time"])
ax1.set_xlabel("Clustering time (s)")
plt.tight_layout()

# %%
# :class:`~sklearn.cluster.KMeans` and :class:`~sklearn.cluster.MiniBatchKMeans`
# suffer from the phenomenon called the `Curse of Dimensionality
# <https://en.wikipedia.org/wiki/Curse_of_dimensionality>`_ for high dimensional
# datasets such as text data. That is the reason why the overall scores improve
# when using LSA. The Silhouette Coefficient is particularly low because its
# definition requires measuring distances, in contrast with other measures such
# as the V-measure and the Adjusted Rand Index which are only based on cluster
# assignments rather than distances.
#
# Using the `"k-means++"` init on LSA reduced data is both good in terms of
# metrics and stability as shown in the plot above. It also has low clustering
# time though the LSA step itself takes a long time, especially with hashed
# vectors. The reason is that a hashed space is typically large (set to
# `n_features=50_000` in this example). One can try lowering the number of
# features at the expense of having a larger fraction of features with hash
# collisions as shown in the example notebook
# :ref:`sphx_glr_auto_examples_text_plot_hashing_vs_dict_vectorizer.py`.
#
# :class:`~sklearn.cluster.MiniBatchKMeans` is unstable for this relatively
# small dataset. It is more interesting to use when the number of samples is
# much bigger, but it comes at the expense of clustering quality.
