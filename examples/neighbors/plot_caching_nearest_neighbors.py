"""
=========================
Caching nearest neighbors
=========================

This examples demonstrates how to precompute the k nearest neighbors before
using them in KNeighborsClassifier. KNeighborsClassifier can compute the
nearest neighbors internally, but precomputing them can have several benefits,
such as finer parameter control, caching for multiple use, or custom
implementations.

Here we use the caching property of pipelines to cache the nearest neighbors
graph between multiple fits of KNeighborsClassifier. The first call is slow
since it computes the neighbors graph, while subsequent call are faster as they
do not need to recompute the graph. Here the durations are small since the
dataset is small, but the gain can be more substantial when the dataset grows
larger, or when the grid of parameter to search is large.

Sample output
-------------
n_neighbors = 05, 0.260 sec
n_neighbors = 10, 0.004 sec
n_neighbors = 15, 0.004 sec
n_neighbors = 20, 0.004 sec
-------------
"""
# Author: Tom Dupre la Tour
#
# License: BSD 3 clause
import time

from sklearn.neighbors import KNeighborsTransformer, KNeighborsClassifier
from sklearn.datasets import load_digits
from sklearn.pipeline import Pipeline

print(__doc__)

X, y = load_digits(return_X_y=True)

n_neighbors_list = [5, 10, 15, 20]

# The graph model uses the maximum number of neighbors necessary in the grid
# search, and the classifier model uses only the nearestas necessary.
graph_model = KNeighborsTransformer(n_neighbors=max(n_neighbors_list),
                                    mode='distance')
classifier_model = KNeighborsClassifier(metric='precomputed')

# Note that we give `memory` a directory to cache the graph computation.
full_model = Pipeline(
    steps=[('graph', graph_model), ('classifier', classifier_model)],
    memory='./cache')

for n_neighbors in n_neighbors_list:
    full_model.named_steps['classifier'].set_params(n_neighbors=n_neighbors)

    start = time.time()
    full_model.fit(X, y)
    print('n_neighbors = %02d, %.3f sec' % (n_neighbors, time.time() - start))
