"""
================================================================
Biclustering documents with the Spectral Co-clustering algorithm
================================================================

This example demonstrates the Spectral Co-clustering algorithm on
the twenty newsgroups dataset.

The TF-IDF vectorized posts form a word frequency matrix, which is
then biclustered using Dhillon's Spectral Co-Clustering algorithm. The
resulting document-word biclusters indicate subsets words used more
often in those subsets documents.

For each bicluster, the category percentage of its documents and its
ten most common words get printed.

The document clusters achieve a better V-measure score than clusters
found by MiniBatchKMeans.

"""
from __future__ import print_function

print(__doc__)

from time import time
import re
from collections import Counter

import numpy as np

from sklearn.datasets.twenty_newsgroups import fetch_20newsgroups
from sklearn.feature_extraction.text import TfidfVectorizer
from sklearn.cluster.bicluster import SpectralCoclustering
from sklearn.cluster import MiniBatchKMeans
from sklearn.metrics.cluster import v_measure_score


def number_aware_tokenizer(doc):
    token_pattern = re.compile(u'(?u)\\b\\w\\w+\\b')
    tokens = token_pattern.findall(doc)
    tokens = ["#NUMBER" if token[0] in "0123456789_" else token
              for token in tokens]
    return tokens


vectorizer = TfidfVectorizer(stop_words='english', min_df=5,
                             tokenizer=number_aware_tokenizer)
cocluster = SpectralCoclustering(n_clusters=20, svd_method='arpack',
                                 random_state=0)
kmeans = MiniBatchKMeans(n_clusters=20, batch_size=5000, random_state=0)

twenty = fetch_20newsgroups()
y_true = twenty.target

print("Vectorizing...")
X = vectorizer.fit_transform(twenty.data)

print("Coclustering...")
start_time = time()
cocluster.fit(X)
y_cocluster = cocluster.row_labels_
print("Done in {}s. V-measure: {}".format(
    time() - start_time,
    v_measure_score(y_cocluster, y_true)))

print("MiniBatchKMeans...")
start_time = time()
y_kmeans = kmeans.fit_predict(X)
print("Done in {}s. V-measure: {}".format(
    time() - start_time,
    v_measure_score(y_kmeans, y_true)))

feature_names = vectorizer.get_feature_names()
document_names = list(twenty.target_names[i] for i in twenty.target)

print("")
print("Biclusters:")
print("----------")
for cluster in xrange(20):
    n_rows, n_cols = cocluster.get_shape(cluster)
    cluster_docs, cluster_words = cocluster.get_indices(cluster)
    if not len(cluster_docs) or not len(cluster_words):
        continue

    # categories
    categories = list(document_names[i] for i in cluster_docs)
    counter = Counter(categories)
    cat_string = ", ".join("{:.0f}% {}".format(float(c) / n_rows * 100,
                                               name)
                           for name, c in counter.most_common()[:3])

    # words
    out_of_cluster_docs = np.where(cocluster.row_labels_[cluster] != cluster)[0]
    word_col = X[:, cluster_words]
    word_scores = np.array(word_col[cluster_docs, :].sum(axis=0) -
                           word_col[out_of_cluster_docs, :].sum(axis=0))
    word_scores = word_scores.ravel()
    important_words = list(feature_names[cluster_words[i]]
                           for i in word_scores.argsort()[:-10:-1])

    print("bicluster {} : {} documents, {} words".format(
        cluster, n_rows, n_cols))
    print("categories  : {}".format(cat_string))
    print("words       : {}\n".format(', '.join(important_words)))
