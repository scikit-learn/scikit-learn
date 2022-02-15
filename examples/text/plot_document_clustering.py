"""
=======================================
Clustering text documents using k-means
=======================================

This is an example showing how the scikit-learn can be used to cluster
documents by topics using a bag-of-words approach. This example uses
a scipy.sparse matrix to store the features instead of standard numpy arrays.
We attempt reweighting and dimensionality reduction on the
extracted features and compare performance metrics for each case.
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

print("%d documents" % len(dataset.data))
print("%d categories" % len(dataset.target_names))
print()


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

from sklearn.decomposition import TruncatedSVD
from sklearn.feature_extraction.text import TfidfVectorizer
from sklearn.pipeline import make_pipeline
from sklearn.preprocessing import Normalizer
from sklearn import metrics
from sklearn.cluster import MiniBatchKMeans


def kmeans_pipeline(data, labels, n_clusters,
                   idf=True, n_components=False, random_state=None):
    pipeline = []

    vectorizer = TfidfVectorizer(
        max_df=0.5,
        max_features=10000,
        min_df=2,
        stop_words="english",
        use_idf=idf,
    )

    X = vectorizer.fit_transform(data)

    pipeline.append(vectorizer)

    if type(n_components) is int:
        # Vectorizer results are normalized, which makes KMeans behave as
        # spherical k-means for better results. Since LSA/SVD results are
        # not normalized, we have to redo the normalization.
        svd = TruncatedSVD(10, random_state=random_state)
        normalizer = Normalizer(copy=False)
        lsa = make_pipeline(svd, normalizer)

        X = lsa.fit_transform(X)

        pipeline.append(svd)

    km = MiniBatchKMeans(
        n_clusters=n_clusters,
        init="k-means++",
        n_init=1,
        init_size=1000,
        batch_size=1000,
        verbose=False,
        random_state=random_state
    )

    km.fit(X)

    pipeline.append(km)

    homo = metrics.homogeneity_score(labels, km.labels_)
    complt = metrics.completeness_score(labels, km.labels_)
    v_msr = metrics.v_measure_score(labels, km.labels_)
    rand_scr = metrics.adjusted_rand_score(labels, km.labels_)
    silhouette = metrics.silhouette_score(X, km.labels_, sample_size=1000)
    print(
        f"{homo:.3f}\t{complt:.3f}\t{v_msr:.3f}\t{rand_scr:.3f}\t{silhouette:.3f}"
        )

    return pipeline


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

labels = dataset.target
true_k = np.unique(labels).shape[0]

rseed = 16

print("         \thomo\tcomplet\tv-meas\trand-i\tsilhouette")
print("-" * 58)

print(f"{'Tf':9s}", end="\t")
kmeans_pipeline(data=dataset.data, labels=labels, n_clusters=true_k,
               idf=False, random_state=rseed)

print(f"{'Tfidf':9s}", end="\t")
kmeans_pipeline(data=dataset.data, labels=labels, n_clusters=true_k,
               idf=True, random_state=rseed)

print(f"{'Tfidf+LSA':9s}", end="\t")
pipeline = kmeans_pipeline(data=dataset.data, labels=labels, n_clusters=true_k,
                          idf=True, n_components=10, random_state=rseed)


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


vectorizer, svd, km = pipeline

original_space_centroids = svd.inverse_transform(km.cluster_centers_)
order_centroids = original_space_centroids.argsort()[:, ::-1]

terms = vectorizer.get_feature_names_out()

fig, ax = plt.subplots(figsize=(10,5))

for i in range(true_k):
    x_pos = .1 + .266 * i
    ax.text(x_pos, .95, f"Cluster {i}",
            fontsize=15, fontweight='bold',
            ha='center', va='bottom')
    for j, ind in enumerate(order_centroids[i, :10]):
        scale = np.log(original_space_centroids[i, ind]+1)
        y_pos = .77-.1*j
        ax.text(x_pos, y_pos, terms[ind],
                color=plt.cm.Set1(i),
                fontsize=scale*130, fontweight=scale*1300,
                ha='center', va='bottom')

ax.set_title("Top terms per cluster", fontsize=20)
ax.axis(False)
plt.show()
