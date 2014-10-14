"""
============================================
Scalability of Approximate Nearest Neighbors
============================================

This demonstrates the relationship between query time and index size of
LSHForest. Also, query time is compared with brute force method in exact
nearest neighbor search for the same index sizes.

Moreover, precision of first 10 neighbors for different index sizes is
also demontrated, LSHForest is built with 'n_estimators' = 3 and
'n_candidates' = 100.
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

# Initialize the range of `n_samples`
n_samples_values = [1000, 2000, 4000, 6000, 8000, 10000]
n_queries = 100
average_times_approx = []
average_times_exact = []
accuracies = []
n_iter = 30

rng = np.random.RandomState(42)

X_all, _ = make_blobs(n_samples=max(n_samples_values) + n_queries,
                      n_features=100, centers=100, shuffle=True,
                      random_state=0)
# Calculate the average query time
for n_samples in n_samples_values:
    X = X_all[:n_samples]
    # Initialize LSHForest for queries of a single neighbor
    lshf = LSHForest(n_estimators=15, n_candidates=200,
                     n_neighbors=10).fit(X)
    nbrs = NearestNeighbors(algorithm='brute', metric='cosine',
                            n_neighbors=10).fit(X)
    average_time_approx = 0
    average_time_exact = 0
    accuracy = 0

    for i in range(n_iter):
        query = X_all[-n_queries:][rng.randint(0, n_queries)]

        t0 = time.time()
        approx_neighbors = lshf.kneighbors(query,
                                           return_distance=False)
        T = time.time() - t0
        average_time_approx += T

        t0 = time.time()
        exact_neighbors = nbrs.kneighbors(query,
                                          return_distance=False)
        T = time.time() - t0
        average_time_exact += T

        accuracy = accuracy + np.in1d(approx_neighbors,
                                      exact_neighbors).mean()

    average_time_approx /= float(n_iter)
    average_time_exact /= float(n_iter)
    accuracy /= float(n_iter)
    average_times_approx.append(average_time_approx)
    average_times_exact.append(average_time_exact)
    accuracies.append(accuracy)

# Plot average query time against n_samples
plt.figure()
colors = ['r', 'b']
p1 = plt.Rectangle((0, 0), 0.1, 0.1, fc=colors[0])
p2 = plt.Rectangle((0, 0), 0.1, 0.1, fc=colors[1])

labels = ['LSHForest', 'NearestNeighbors']
plt.legend((p1, p2), (labels[0], labels[1]), loc='upper left')

plt.scatter(n_samples_values, average_times_approx, c='r')
plt.plot(n_samples_values, average_times_approx, c='r')
plt.scatter(n_samples_values, average_times_exact, c='b')
plt.plot(n_samples_values, average_times_exact, c='b')
plt.xlim(min(n_samples_values), max(n_samples_values))
plt.ylim(min(average_times_approx + average_times_exact),
         max(average_times_approx + average_times_exact))
plt.ylabel("Average time in seconds")
plt.xlabel("n_samples")
plt.grid(which='both')
plt.title("Impact of index size on response time for first "
          "nearest neighbors queries")

# Plot precision@10
plt.figure()
plt.scatter(n_samples_values, accuracies, c='c')
plt.plot(n_samples_values, accuracies, c='c')
plt.xlim(min(n_samples_values), max(n_samples_values))
plt.ylim(0, 1)
plt.ylabel("precision@10")
plt.xlabel("n_samples")
plt.grid(which='both')
plt.title("precision of first 10 nearest neighbor queries with index size")

plt.show()
