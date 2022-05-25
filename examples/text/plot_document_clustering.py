"""
=======================================
Clustering text documents using k-means
=======================================

This is an example showing how the scikit-learn can be used to cluster
documents by topics using a bag-of-words approach. This example uses
a scipy.sparse matrix to store the features instead of standard numpy arrays.

Two feature extraction methods can be used in this example:

  - TfidfVectorizer uses a in-memory vocabulary (a python dict) to map the most
    frequent words to features indices and hence compute a word occurrence
    frequency (sparse) matrix. The word frequencies are then reweighted using
    the Inverse Document Frequency (IDF) vector collected feature-wise over
    the corpus.

  - HashingVectorizer hashes word occurrences to a fixed dimensional space,
    possibly with collisions. The word count vectors are then normalized to
    each have l2-norm equal to one (projected to the euclidean unit-ball) which
    seems to be important for k-means to work in high dimensional space.

    HashingVectorizer does not provide IDF weighting as this is a stateless
    model (the fit method does nothing). When IDF weighting is needed it can
    be added by pipelining its output to a TfidfTransformer instance.

Two algorithms are demoed: ordinary k-means and its more scalable cousin
minibatch k-means.

Additionally, latent semantic analysis can also be used to reduce
dimensionality and discover latent patterns in the data.

It can be noted that k-means (and minibatch k-means) are very sensitive to
feature scaling and that in this case the IDF weighting helps improve the
quality of the clustering by quite a lot as measured against the "ground truth"
provided by the class label assignments of the 20 newsgroups dataset.

This improvement is not visible in the Silhouette Coefficient which is small
for both as this measure seem to suffer from the phenomenon called
"Concentration of Measure" or "Curse of Dimensionality" for high dimensional
datasets such as text data. Other measures such as V-measure and Adjusted Rand
Index are information theoretic based evaluation scores: as they are only based
on cluster assignments rather than distances, hence not affected by the curse
of dimensionality.

Note: as k-means is optimizing a non-convex objective function, it will likely
end up in a local optimum. Several runs with independent random init might be
necessary to get a good convergence.

"""

# Author: Peter Prettenhofer <peter.prettenhofer@gmail.com>
#         Lars Buitinck
# License: BSD 3 clause


# %%
# Load Data
# =========
#
# We load data from the :ref:`20newsgroups_dataset`, which comprises around
# 18000 newsgroups posts on 20 topics. For illustrative purposes and reducing
# the computational cost, we select a subset of 4 topics only.  See the example
# :ref:`sphx_glr_auto_examples_text_plot_document_classification_20newsgroups.py`
# to gain intuition on the overlap of such topics.

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
# Feature Extraction using sparse vectorizers
# ===========================================

from sklearn.feature_extraction.text import HashingVectorizer
from sklearn.feature_extraction.text import TfidfTransformer
from sklearn.pipeline import make_pipeline
from time import time

# Perform an IDF normalization on the output of HashingVectorizer

hasher = HashingVectorizer(
    stop_words="english",
    alternate_sign=False,
    norm=None,
)
vectorizer = make_pipeline(hasher, TfidfTransformer())

t0 = time()
X = vectorizer.fit_transform(dataset.data)
print(f"done in {time() - t0:.3f}s")
print(f"n_samples: {X.shape[0]}, n_features: {X.shape[1]}")

# %%
# Performing dimensionality reduction using LSA
# Since LSA/SVD results are not normalized,
# we redo the normalization to improve the k-means result.

from sklearn.decomposition import TruncatedSVD
from sklearn.preprocessing import Normalizer

t0 = time()

svd = TruncatedSVD(n_components=100)
normalizer = Normalizer(copy=False)
lsa = make_pipeline(svd, normalizer)

X = lsa.fit_transform(X)

print(f"done in {time() - t0:.3f}s")

explained_variance = svd.explained_variance_ratio_.sum()
print(f"Explained variance of the SVD step: {int(explained_variance * 100)}%")

print()

# %%
# TfidfVectorizer
# ---------------
from sklearn.feature_extraction.text import TfidfVectorizer

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
# Performing dimensionality reduction using LSA
# Since LSA/SVD results are not normalized,
# we redo the normalization to improve the k-means result.

t0 = time()
svd = TruncatedSVD(n_components=100)
normalizer = Normalizer(copy=False)
lsa = make_pipeline(svd, normalizer)
X = lsa.fit_transform(X)

print(f"done in {time() - t0}s")

explained_variance = svd.explained_variance_ratio_.sum()
print("Explained variance of the SVD step: {}%".format(int(explained_variance * 100)))

# %%
# Clustering
# ==========

from sklearn.cluster import KMeans, MiniBatchKMeans


km = MiniBatchKMeans(
    n_clusters=true_k,
    init="k-means++",
    n_init=1,
    init_size=1000,
    batch_size=1000,
    # verbose=opts.verbose,
)

print("Clustering sparse data with %s" % km)
t0 = time()
km.fit(X)
print("done in %0.3fs" % (time() - t0))
print()

km = KMeans(
    n_clusters=true_k,
    init="k-means++",
    max_iter=100,
    n_init=1,
    # verbose=opts.verbose,
)

print("Clustering sparse data with %s" % km)
t0 = time()
km.fit(X)
print("done in %0.3fs" % (time() - t0))
print()


# %%
# Performance metrics
# ===================

from sklearn import metrics

print("Homogeneity: %0.3f" % metrics.homogeneity_score(labels, km.labels_))
print("Completeness: %0.3f" % metrics.completeness_score(labels, km.labels_))
print("V-measure: %0.3f" % metrics.v_measure_score(labels, km.labels_))
print("Adjusted Rand-Index: %.3f" % metrics.adjusted_rand_score(labels, km.labels_))
print(
    "Silhouette Coefficient: %0.3f"
    % metrics.silhouette_score(X, km.labels_, sample_size=1000)
)

print()


# %%
# Top terms per cluster


original_space_centroids = svd.inverse_transform(km.cluster_centers_)
order_centroids = original_space_centroids.argsort()[:, ::-1]


terms = vectorizer.get_feature_names_out()
for i in range(true_k):
    print("Cluster %d:" % i, end="")
    for ind in order_centroids[i, :10]:
        print(" %s" % terms[ind], end="")
    print()

# %%
