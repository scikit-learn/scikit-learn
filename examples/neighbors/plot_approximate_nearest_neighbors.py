"""
=============================
Approximate Nearest Neighbors
=============================

Locality Sensitive Hashing is used as the approximate nearest neighbor
search method in this example. This demonstrates the behaviour of the
accuracy of the nearest neighbor queries of Locality Sensitive Hashing
Forest as the number of candidates and the number of estimators (trees)
vary.

In the first plot, the number of trees is fixed (kept at default value 10).
Accuracy is measured with the number of candidates. Here, the term
"number of candidates" refers to maximum bound for the number of distinct
points retrieved from each tree to calculate the distances. Nearest
neighbors are selected from this pool of candidates.

In the second plot, the number of candidates is fixed at 500. Number of trees
is varied and the accuracy is plotted against those values. To measure the
accuracy, the true nearest neighbors are required, therefore
:class:`sklearn.neighbors.NearestNeighbors` is used to compute the exact
neighbors.

Third plot demonstrates the behavior of approximate nearest neighbor queries
of LSHForest as the size of the dataset varies. Here, average query time for
a single neighbor vs number of samples is plotted.
"""
from __future__ import division
print(__doc__)

# Author: Maheshakya Wijewardena <maheshakya.10@cse.mrt.ac.lk>
#
# License: BSD 3 clause


###############################################################################
import time
import numpy as np
from sklearn.datasets.samples_generator import make_blobs
from sklearn.neighbors import LSHForest
from sklearn.neighbors import NearestNeighbors
import matplotlib.pyplot as plt


# Initialize size of the database, iterations and required neighbors.
n_samples = 10000
n_features = 100
n_queries = 30
rng = np.random.RandomState(42)

# Generate sample data
X, _ = make_blobs(n_samples=n_samples + n_queries,
                  n_features=n_features, centers=10,
                  cluster_std=5, random_state=0)
X_index = X[:n_samples]
X_query = X[n_samples:]
# Set `n_candidate` values
n_candidates_values = np.linspace(10, 500, 5).astype(np.int)
n_estimators_for_candidate_value = [1, 5, 10]
accuracies_c = np.zeros((len(n_estimators_for_candidate_value),
                         n_candidates_values.shape[0]), dtype=float)

# Calculate average accuracy for each value of `n_candidates`
for j, value in enumerate(n_estimators_for_candidate_value):
    for i, n_candidates in enumerate(n_candidates_values):
        lshf = LSHForest(n_estimators=value,
                         n_candidates=n_candidates, n_neighbors=1)
        nbrs = NearestNeighbors(n_neighbors=1)
        # Fit the Nearest neighbor models
        lshf.fit(X_index)
        nbrs.fit(X_index)
        # Get neighbors
        neighbors_approx = lshf.kneighbors(X_query, return_distance=False)
        neighbors_exact = nbrs.kneighbors(X_query, return_distance=False)

        accuracies_c[j, i] = np.sum(np.equal(neighbors_approx,
                                             neighbors_exact))/n_queries

# Set `n_estimators` values
n_estimators_values = np.arange(1, 11, 1).astype(np.int)
accuracies_trees = np.zeros(n_estimators_values.shape[0], dtype=float)

# Calculate average accuracy for each value of `n_estimators`
for i, n_estimators in enumerate(n_estimators_values):
    lshf = LSHForest(n_estimators=n_estimators, n_neighbors=1)
    nbrs = NearestNeighbors(n_neighbors=1)
    # Fit the Nearest neighbor models
    lshf.fit(X_index)
    nbrs.fit(X_index)
    # Get neighbors
    neighbors_approx = lshf.kneighbors(X_query, return_distance=False)
    neighbors_exact = nbrs.kneighbors(X_query, return_distance=False)

    accuracies_trees[i] = np.sum(np.equal(neighbors_approx,
                                          neighbors_exact))/n_queries

###############################################################################
# Plot the accuracy variation with `n_estimators`
plt.figure()
plt.scatter(n_estimators_values, accuracies_trees, c='k')
plt.plot(n_estimators_values, accuracies_trees, c='g')
plt.ylim([0, 1])
plt.xlim(min(n_estimators_values), max(n_estimators_values))
plt.ylabel("Accuracy")
plt.xlabel("n_estimators")
plt.grid(which='both')
plt.title("Accuracy variation with n_estimators")

plt.show()

# Plot the accuracy variation with `n_candidates`
colors = ['c', 'm', 'y']
p1 = plt.Rectangle((0, 0), 0.1, 0.1, fc=colors[0])
p2 = plt.Rectangle((0, 0), 0.1, 0.1, fc=colors[1])
p3 = plt.Rectangle((0, 0), 0.1, 0.1, fc=colors[2])

labels = ['n_estimators = ' + str(n_estimators_for_candidate_value[0]),
          'n_estimators = ' + str(n_estimators_for_candidate_value[1]),
          'n_estimators = ' + str(n_estimators_for_candidate_value[2])]

plt.figure()
plt.legend((p1, p2, p3), (labels[0], labels[1], labels[2]), loc='upper left')

for i in range(len(n_estimators_for_candidate_value)):
    plt.scatter(n_candidates_values, accuracies_c[i, :], c=colors[i])
    plt.plot(n_candidates_values, accuracies_c[i, :], c=colors[i])
plt.ylim([0, 1])
plt.xlim(min(n_candidates_values), max(n_candidates_values))
plt.ylabel("Accuracy")
plt.xlabel("n_candidates")
plt.grid(which='both')
plt.title("Accuracy variation with n_candidates")

plt.show()

###############################################################################
# Initialize the range of `n_samples`
n_samples_values = [10, 100, 1000, 10000, 100000]
average_times = []
n_iter = 30
# Calculate the average query time
for n_samples in n_samples_values:
    X, labels_true = make_blobs(n_samples=n_samples, n_features=10,
                                centers=10, cluster_std=5,
                                random_state=0)
    # Initialize LSHForest for queries of a single neighbor
    lshf = LSHForest(n_candidates=1000, n_neighbors=1)
    lshf.fit(X)

    average_time = 0

    for i in range(n_iter):
        query = X[rng.randint(0, n_samples)]
        t0 = time.time()
        approx_neighbors = lshf.kneighbors(query,
                                           return_distance=False)[0]
        T = time.time() - t0
        average_time = average_time + T

    average_time = average_time / float(n_iter)
    average_times.append(average_time)

# Plot average query time against n_samples
plt.scatter(n_samples_values, average_times, c='k')
plt.plot(n_samples_values, average_times, c='b')
plt.xlim(min(n_samples_values), max(n_samples_values))
plt.ylim(min(average_times), max(average_times))
plt.semilogx()
plt.ylabel("Average time in seconds")
plt.xlabel("n_samples")
plt.title("Average time vs n_samples")

plt.show()
