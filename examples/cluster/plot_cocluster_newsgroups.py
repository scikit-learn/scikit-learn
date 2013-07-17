"""
====================================
Spectral Co-clustering of newsgroups
====================================

This example demonstrates the spectral co-clustering algorithm on a
few categories from the twenty newsgroups dataset.

The TF-IDF vectorized posts form a word frequency matrix, which is
then biclustered using Dhillon's Spectral Co-Clustering algorithm. The
resulting document-word biclusters indicate subsets words used more
often in those subsets documents.

Another way to look at this result: given a document clustering, the
words could be clustered to associate with each document cluster the
words used more often in that cluster than in any other. Similarly,
given a word clustering, the documents could be clustered in the same
way. The Spectral Co-Clustering algorithm does both simultaneously.

For each bicluster, the category percentage of its documents and its
ten most common words get printed.

"""
print(__doc__)

# Author: Kemal Eren <kemal@kemaleren.com>
# License: BSD 3 clause

from collections import Counter

import numpy as np

from sklearn.datasets import fetch_20newsgroups
from sklearn.feature_extraction.text import TfidfVectorizer
from sklearn.cluster.bicluster import SpectralCoclustering

categories = ['alt.atheism', 'comp.graphics']
newsgroups = fetch_20newsgroups(subset='all', shuffle=False,
                                categories=categories)

vec = TfidfVectorizer(stop_words='english', min_df=0.01)
data = vec.fit_transform(newsgroups.data)

row_names = list(newsgroups.target_names[i] for i in newsgroups.target)
feature_names = vec.get_feature_names()

model = SpectralCoclustering(n_clusters=len(categories),
                             svd_method='randomized', random_state=0)
model.fit(data)

for i in range(model.rows_.shape[0]):
    n_rows, n_cols = model.get_shape(i)
    row_idx, col_idx = model.get_indices(i)

    # categories
    category_indices = list(row_names[i] for i in row_idx)
    counter = Counter(category_indices)

    cat_string = ", ".join("{:.0f}% {}".format(float(c) / n_rows * 100,
                                               name)
                           for name, c in counter.most_common())

    # words
    word_sum = data[row_idx, :][:, col_idx].sum(axis=0)
    word_sum = np.asarray(word_sum).squeeze()
    important_indices = np.argsort(word_sum)[::-1][:10]
    words = list(feature_names[col_idx[i]] for i in important_indices)

    print "bicluster {} : {} documents, {} words".format(i, n_rows, n_cols)
    print "categories  : {}".format(cat_string)
    print "words       : {}\n".format(', '.join(words))
