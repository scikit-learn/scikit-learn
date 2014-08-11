"""
=================================
Locality Sensitive Hashing Forest
=================================

Demonstrate the behaviour of the accuracy of the nearest neighbor
queries of Locality Sensitive Hashing Forest as the number of candidates
and the number of estimators(trees) vary.

"""
print(__doc__)

# Author: Maheshakya Wijewardena <maheshakya.10@cse.mrt.ac.lk>
#
# License: BSD 3 clause


###############################################################################
import numpy as np
from sklearn.metrics import euclidean_distances
from sklearn.datasets.samples_generator import make_blobs
from sklearn.neighbors import LSHForest

# Initialize size of the database, iterations and required neighbors.
n_samples = 10000
n_dim = 100
n_iter = 100
n_neighbors = 100

# Generate sample data
X, _ = make_blobs(n_samples=n_samples, n_features=n_dim,
                  centers=10, cluster_std=5, random_state=0)

# Set `n_candidate` values
n_candidates_values = np.arange(10, 1010, 100)
accuracies_c = np.zeros(n_candidates_values.shape[0], dtype=float)

# Calculate average accuracy for each value of `n_candidates`
for i, n_candidates in enumerate(n_candidates_values):
    lshf = LSHForest(n_candidates=n_candidates, n_neighbors=n_neighbors)
    # Fit the LSHForest model
    lshf.fit(X)
    for j in range(n_iter):
        point = X[np.random.randint(0,n_samples)]
        # Get neighbors
        neighbors = lshf.kneighbors(point)
        distances = euclidean_distances(point, X)
        ranks = np.argsort(distances)[0,:n_neighbors]

        intersection = np.intersect1d(ranks, neighbors).shape[0]
        ratio = intersection/float(n_neighbors)
        accuracies_c[i] = accuracies_c[i] + ratio

    accuracies_c[i] = accuracies_c[i]/float(n_iter)

# Set `n_estimators` values
n_estimators_values = np.arange(1, 100, 5)
accuracies_trees = np.zeros(n_estimators_values.shape[0], dtype=float)

# Calculate average accuracy for each value of `n_estimators`
for i, n_estimators in enumerate(n_estimators_values):
    lshf = LSHForest(n_candidates=500, n_estimators=n_estimators, n_neighbors=n_neighbors)
    lshf.fit(X)
    for j in range(n_iter):
        point = X[np.random.randint(0,n_samples)]
        neighbors = lshf.kneighbors(point)
        distances = euclidean_distances(point, X)
        ranks = np.argsort(distances)[0,:n_neighbors]

        intersection = np.intersect1d(ranks, neighbors).shape[0]
        ratio = intersection/float(n_neighbors)
        accuracies_trees[i] = accuracies_trees[i] + ratio

    accuracies_trees[i] = accuracies_trees[i]/float(n_iter)

###############################################################################
# Plot the accuracy resutls
plt.subplot(2, 1, 0)
plt.scatter(n_candidates_values, accuracies_c, c='k')
plt.plot(n_candidates_values, accuracies_c, c='g')
plt.ylim([0,1.1])
plt.xlim([0,1000])
plt.ylabel("Accuracy")
plt.legend()
plt.title("Accuracy variation with n_candidates")

plt.subplot(2, 1, 1)
plt.scatter(n_estimators_values, accuracies_trees, c='k')
plt.plot(n_estimators_values, accuracies_trees, c='r')
plt.ylim([0,1.1])
plt.xlim([0,100])
plt.ylabel("Accuracy")
plt.legend()
plt.title("Accuracy variation with n_estimators")

plt.show()