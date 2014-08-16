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
`sklearn.neighbors.NearestNeighbors` is used to compute the exact neighbors.

Third plot demonstrates the behavior of approximate nearest neighbor queries
of LSHForest as the size of the dataset varies. Here, average query time for
a single neighbor vs number of samples is plotted.
"""
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
n_iter = 30
n_neighbors = 100
rng = np.random.RandomState(42)

# Generate sample data
X, _ = make_blobs(n_samples=n_samples, n_features=n_features,
                  centers=10, cluster_std=5, random_state=0)

# Set `n_candidate` values
n_candidates_values = np.linspace(10, 500, 5).astype(np.int)
accuracies_c = np.zeros(n_candidates_values.shape[0], dtype=float)

# Calculate average accuracy for each value of `n_candidates`
for i, n_candidates in enumerate(n_candidates_values):
    lshf = LSHForest(n_candidates=n_candidates, n_neighbors=n_neighbors)
    nbrs = NearestNeighbors(n_neighbors=n_neighbors, algorithm='brute')
    # Fit the Nearest neighbor models
    lshf.fit(X)
    nbrs.fit(X)
    for j in range(n_iter):
        query = X[rng.randint(0, n_samples)]
        # Get neighbors
        neighbors_approx = lshf.kneighbors(query, return_distance=False)
        neighbors_exact = nbrs.kneighbors(query, return_distance=False)

        intersection = np.intersect1d(neighbors_approx,
                                      neighbors_exact).shape[0]
        ratio = intersection / float(n_neighbors)
        accuracies_c[i] += ratio

    accuracies_c[i] = accuracies_c[i] / float(n_iter)

# Set `n_estimators` values
n_estimators_values = np.linspace(1, 30, 5).astype(np.int)
accuracies_trees = np.zeros(n_estimators_values.shape[0], dtype=float)

# Calculate average accuracy for each value of `n_estimators`
for i, n_estimators in enumerate(n_estimators_values):
    lshf = LSHForest(n_candidates=500, n_estimators=n_estimators,
                     n_neighbors=n_neighbors)
    nbrs = NearestNeighbors(n_neighbors=n_neighbors, algorithm='brute')

    lshf.fit(X)
    nbrs.fit(X)
    for j in range(n_iter):
        query = X[rng.randint(0, n_samples)]
        neighbors_approx = lshf.kneighbors(query, return_distance=False)
        neighbors_exact = nbrs.kneighbors(query, return_distance=False)

        intersection = np.intersect1d(neighbors_approx,
                                      neighbors_exact).shape[0]
        ratio = intersection / float(n_neighbors)
        accuracies_trees[i] += ratio

    accuracies_trees[i] = accuracies_trees[i] / float(n_iter)

###############################################################################
# Plot the accuracy variation with `n_estimators`
plt.subplot(2, 1, 0)
plt.scatter(n_candidates_values, accuracies_c, c='k')
plt.plot(n_candidates_values, accuracies_c, c='g')
plt.ylim([0, 1])
plt.xlim(min(n_candidates_values), max(n_candidates_values))
plt.ylabel("Accuracy")
plt.title("Accuracy variation with n_candidates")

# Plot the accuracy variation with `n_candidates`
plt.subplot(2, 1, 1)
plt.scatter(n_estimators_values, accuracies_trees, c='k')
plt.plot(n_estimators_values, accuracies_trees, c='r')
plt.ylim([0, 1])
plt.xlim(min(n_estimators_values), max(n_estimators_values))
plt.ylabel("Accuracy")
plt.title("Accuracy variation with n_estimators")

plt.show()

###############################################################################
# Initialize the range of `n_samples`
n_samples_values = [10, 100, 1000, 10000, 100000]
average_times = []
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
