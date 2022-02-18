"""
=======================================
Clustering text documents using k-means
=======================================

This is an example showing how scikit-learn can be used to cluster
documents by topics using a bag-of-words approach. This example uses
a `scipy.sparse` matrix to store the features instead of standard numpy arrays.
Here, we compare different feature extraction approaches using several
clustering metrics.
"""

# Author: Peter Prettenhofer <peter.prettenhofer@gmail.com>
#         Lars Buitinck
#         Giorgos Papadokostakis
# License: BSD 3 clause


# %%
# Load some categories from the training set
# ------------------------------------------
# We will start by loading 4 of the 20 news groups categories and see if
# k-means can cluster them properly.

from sklearn.datasets import fetch_20newsgroups
import numpy as np


categories = [
    "alt.atheism",
    "talk.religion.misc",
    "comp.graphics",
    "sci.space",
]

print("Loading 20 newsgroups dataset for categories:")
print(categories)

dataset = fetch_20newsgroups(
    subset="all", categories=categories, shuffle=True, random_state=42
)

data = dataset.data
labels = dataset.target
true_k = np.unique(labels).shape[0]

print(f"{len(data)} documents")
print(f"{len(dataset.target_names)} categories")


# %%
# Define k-means pipeline
# -----------------------
# Here we define our k-means pipeline and the metrics to evaluate performance.
#
# Term-frequency inverse document-frequency (Tfidf) is used as a feature
# extraction method.
# TfidfVectorizer uses an in-memory vocabulary (a python dict) to map the most
# frequent words to features indices and hence compute a word occurrence
# frequency (sparse) matrix. The word frequencies are then reweighted using
# the Inverse Document Frequency (IDF) vector collected feature-wise over
# the corpus.
#
# Additionally, latent semantic analysis (LSA) is also be used to reduce
# dimensionality and discover latent patterns in the data.
#
# Instead of ordinary k-means we use its more scalable cousin
# minibatch k-means.


# %%

from sklearn.feature_extraction.text import TfidfVectorizer
from sklearn.cluster import KMeans


tf_model = [
    TfidfVectorizer(
        max_df=0.5,
        max_features=10000,
        min_df=2,
        stop_words="english",
        use_idf=False,
    ),
    KMeans(n_clusters=true_k, init="k-means++", random_state=0),
]


# %%

tf_idf_model = [
    TfidfVectorizer(
        max_df=0.5,
        max_features=10000,
        min_df=2,
        stop_words="english",
        use_idf=True,
    ),
    KMeans(n_clusters=true_k, init="k-means++", random_state=0),
]


# %%

from sklearn.decomposition import TruncatedSVD
from sklearn.preprocessing import Normalizer


tf_idf_lsa_model = [
    TfidfVectorizer(
        max_df=0.5,
        max_features=10000,
        min_df=2,
        stop_words="english",
        use_idf=True,
    ),
    TruncatedSVD(n_components=10, random_state=0),
    Normalizer(copy=False),
    KMeans(n_clusters=true_k, init="k-means++", random_state=0),
]


# %%

from collections import defaultdict

import pandas as pd

from sklearn.pipeline import make_pipeline
from sklearn.utils import resample
from sklearn.metrics import (
    homogeneity_score,
    completeness_score,
    v_measure_score,
    adjusted_rand_score,
    silhouette_score,
)


def fit_evaluate_clusterer_bootstrap(
    model, X, labels, n_bootstrap=10, random_state=None
):
    bootstrap_results = defaultdict(list)
    rng = np.random.RandomState(random_state)
    vectorizer, km = make_pipeline(*model[:-1]), model[-1]

    for _ in range(n_bootstrap):
        X_bootstrap, labels_bootstrap = resample(
            X, labels, replace=True, random_state=rng
        )

        X_vec = vectorizer.fit_transform(X_bootstrap)
        predicted_labels = km.fit_predict(X_vec)

        bootstrap_results["homogeneity"].append(
            homogeneity_score(labels_bootstrap, predicted_labels)
        )
        bootstrap_results["completeness"].append(
            completeness_score(labels_bootstrap, predicted_labels)
        )
        bootstrap_results["v-measure"].append(
            v_measure_score(labels_bootstrap, predicted_labels)
        )
        bootstrap_results["ARI"].append(
            adjusted_rand_score(labels_bootstrap, predicted_labels)
        )
        bootstrap_results["silhouette"].append(
            silhouette_score(X_vec, predicted_labels, sample_size=20)
        )

    return pd.DataFrame(
        bootstrap_results, index=[f"Bootstrap #{i}" for i in range(n_bootstrap)]
    )


# %%
# Compare performance metrics
# ---------------------------
# It can be noted that k-means (and minibatch k-means) are very sensitive to
# feature scaling and that in this case the IDF weighting helps improve the
# quality of the clustering by quite a lot as measured against the "ground truth"
# provided by the class label assignments of the 20 newsgroups dataset.
# We can see this improvement in metrics homogeneity, completeness, V-measure
# and adjusted Rand index.
#
# This improvement is not visible in the Silhouette Coefficient which is small
# for both as this measure seem to suffer from the phenomenon called
# "Concentration of Measure" or "Curse of Dimensionality" for high dimensional
# datasets such as text data. The other measures are information theoretic
# based evaluation scores: as they are only based
# on cluster assignments rather than distances, hence not affected by the curse
# of dimensionality.
# We can see a significant improvemend though in Silhouette Coefficient with
# LSA as the dimensionality reduction makes the space more compact and
# the Euclidian distances between samples much lower.
#
# Note: as k-means is optimizing a non-convex objective function, it will
# likely end up in a local optimum. Several runs with independent random init
# might be necessary to get a good convergence. Here we have chosen a random
# state that produces good results.

pd.set_option("display.max_columns", 10)

results = []

results.append(
    fit_evaluate_clusterer_bootstrap(
        tf_model, data, labels, n_bootstrap=2, random_state=0
    )
)
results.append(
    fit_evaluate_clusterer_bootstrap(
        tf_idf_model, data, labels, n_bootstrap=2, random_state=0
    )
)
results.append(
    fit_evaluate_clusterer_bootstrap(
        tf_idf_lsa_model, data, labels, n_bootstrap=2, random_state=0
    )
)

results = pd.concat(
    [df.aggregate(["mean", "std"]).melt().transpose() for df in results], axis=0
)

index = pd.MultiIndex.from_product(
    [
        ["homogeneity", "completeness", "v-measure", "ARI", "silhouette"],
        ["mean", "std"],
    ]
)
performance_summary = pd.DataFrame(
    results.loc["value"].values, index=["Tf", "Tf-idf", "Tf-idf-LSA"], columns=index
)
print(performance_summary.round(3))


# %%
# Plot top terms per cluster
# --------------------------
# Finally we attempt to plot the top ten terms for each cluster centroid, with
# their font size scaled by their weighted term frequency.
# Notice how k-means seems to have clustered together two of the four categories
# of the data and created another unrelated cluster.
# This seems to be a common theme of k-means, when the data categories are
# closely related it lumps them together.

import matplotlib.pyplot as plt


vectorizer, svd, _, km = tf_idf_lsa_model

original_space_centroids = svd.inverse_transform(km.cluster_centers_)
order_centroids = original_space_centroids.argsort()[:, ::-1]

terms = vectorizer.get_feature_names_out()

fig, ax = plt.subplots(figsize=(10, 5))

for i in range(true_k):
    x_pos = 0.1 + 0.266 * i
    ax.text(
        x_pos,
        0.95,
        f"Cluster {i}",
        fontsize=15,
        fontweight="bold",
        ha="center",
        va="bottom",
    )
    for j, ind in enumerate(order_centroids[i, :10]):
        scale = np.log(original_space_centroids[i, ind] + 1)
        y_pos = 0.77 - 0.1 * j
        ax.text(
            x_pos,
            y_pos,
            terms[ind],
            color=plt.cm.Set1(i),
            fontsize=scale * 130,
            fontweight=scale * 1300,
            ha="center",
            va="bottom",
        )

ax.set_title("Top terms per cluster", fontsize=20)
ax.axis("off")
plt.show()
