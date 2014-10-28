"""
============================
Nearest Neighbors regression
============================

Demonstrate the resolution of a regression problem
using a k-Nearest Neighbor and the interpolation of the
target using both barycenter and constant weights.

"""
print(__doc__)

# Author: Alexandre Gramfort <alexandre.gramfort@inria.fr>
#         Fabian Pedregosa <fabian.pedregosa@inria.fr>
#
# License: BSD 3 clause (C) INRIA


###############################################################################
# Generate sample data
import numpy as np
import matplotlib.pyplot as plt
from sklearn import neighbors

np.random.seed(0)
X_train = np.sort(5 * np.random.rand(40, 1), axis=0)
X_test = np.linspace(0, 5, 500)[:, np.newaxis]
y_train = np.sin(X_train).ravel()
y_test = np.sin(X_test).ravel()

# Add noise to targets
y_train += 0.2 * np.random.randn(y_train.size)

###############################################################################
# Fit regression model
n_neighbors = 5

plt.figure(figsize=(8, 9))
for i, weights in enumerate(['uniform', 'distance', 'epanechnikov']):
    knn = neighbors.KNeighborsRegressor(n_neighbors, weights=weights)
    y_pred = knn.fit(X_train, y_train).predict(X_test)
    mse = mean_squared_error(y_test, y_pred)
    plt.subplot(3, 1, i + 1)
    plt.scatter(X_train, y_train, c='k', label='training data')
    plt.plot(X_test, y_pred, c='g', label='prediction')
    plt.plot(X_test, y_test, c='k', label='true function')    
    plt.axis('tight')
    plt.legend()
    plt.title("KNeighborsRegressor (k = %i, weights = '%s')" % (n_neighbors,
                                                                weights))    
plt.tight_layout(pad=0.5)
plt.show()
