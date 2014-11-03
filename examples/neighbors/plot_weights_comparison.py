"""
========================================================================
Nearest neighbors classification with different weights and n_neighbors
========================================================================

This example demonstrates how accuracy of k-nearest-neighbors classification
depends on number of neighbors and weights assigned to them. The real-world
Landsat dataset is used. The train and test sets contain 3000 and 3435 samples
respectively.
"""


# Author: Nikolay Mayorov <n59_ru@hotmail.com>
# License: BSD 3 clause


import numpy as np
import matplotlib.pyplot as plt
from sklearn.datasets import fetch_landsat
from sklearn.cross_validation import train_test_split
from sklearn.neighbors import KNeighborsClassifier

landsat = fetch_landsat()
X = landsat.data
y = landsat.target
X_train, X_test, y_train, y_test = train_test_split(X, y, train_size=4000,
                                                    random_state=0)

weight_types = ['uniform', 'distance', 'linear']
results = {w: [] for w in weight_types}
neighbor_counts = np.arange(1, 16, 2)
knn = KNeighborsClassifier()
knn.fit(X_train, y_train)
for weights in weight_types:
    knn.set_params(weights=weights)
    for n_neighbors in neighbor_counts:
        knn.set_params(n_neighbors=n_neighbors)
        results[weights].append(knn.score(X_test, y_test))

for weights in weight_types:
    plt.plot(neighbor_counts, results[weights],
             'o-', label="weights='{}'".format(weights))

plt.xticks(neighbor_counts)
plt.title("Accuracy of kNN on Landsat dataset")
plt.xlabel("n_neighbors")
plt.ylabel("Accuracy")
plt.legend(loc='best')
plt.show()
