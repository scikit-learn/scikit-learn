"""
==========================================
Large Margin Nearest Neighbor Illustration
==========================================

An example illustrating the goal of learning a distance metric that maximizes
the nearest neighbors classification accuracy. The example is solely for
illustration purposes. Please refer to the :ref:`User Guide <lmnn>` for
more information.
"""

# Author: John Chiotellis <johnyc.code@gmail.com>
# License: BSD 3 clause

import numpy as np
import matplotlib.pyplot as plt
from matplotlib.colors import ListedColormap

from sklearn.datasets import make_classification
from sklearn.neighbors import KNeighborsClassifier, \
    NeighborhoodComponentsAnalysis

from sklearn.pipeline import Pipeline

print(__doc__)

n_neighbors = 1
random_state = 0

cmap_light = ListedColormap(['#FFAAAA', '#AAFFAA', '#AAAAFF'])
cmap_bold = ListedColormap(['#FF0000', '#00FF00', '#0000FF'])

# Create a tiny data set of 9 samples from 3 classes
X, y = make_classification(n_samples=9, n_features=2, n_informative=2,
                           n_redundant=0, n_classes=3, n_clusters_per_class=1,
                           class_sep=1.0, random_state=random_state)

X = np.vstack([np.hstack([np.linspace(0, 1, 5)[:, np.newaxis],
                        np.linspace(0, 1, 5)[:, np.newaxis] + 0.3]),
              np.hstack([np.linspace(0, 1, 3)[:, np.newaxis],
                        np.linspace(0, 1, 3)[:, np.newaxis]])])
y = np.hstack([np.zeros(5), np.ones(3)])

n_neighbors = 1

x_min, x_max = X[:, 0].min() - 0.2, X[:, 0].max() + 0.2
y_min, y_max = X[:, 1].min() - 0.2, X[:, 1].max() + 0.2
xx, yy = np.meshgrid(np.arange(x_min, x_max, 0.005),
                     np.arange(y_min, y_max, 0.005))
grid = np.c_[xx.ravel(), yy.ravel()]

# plot decision boundary with knn in the input space
knn = KNeighborsClassifier(n_neighbors)
knn.fit(X, y)

plt.figure()
Z = knn.predict(grid).reshape(xx.shape)
plt.pcolormesh(xx, yy, Z, cmap=cmap_light, alpha=.8)
plt.scatter(X[:, 0], X[:, 1], c=y, cmap=cmap_bold, s=20, edgecolor='k')

# plot decision boundary with NCA in the input space
nca = NeighborhoodComponentsAnalysis(verbose=1)
X_embedded = nca.fit_transform(X, y)

plt.figure()
grid_e = nca.transform(grid)
knn.fit(X_embedded, y)
Z = knn.predict(grid_e).reshape(xx.shape)
plt.pcolormesh(xx, yy, Z, cmap=cmap_light, alpha=.8)
plt.scatter(X[:, 0], X[:, 1], c=y, cmap=cmap_bold, s=20, edgecolor='k')

# # plot decision boundary with nca in the embedding space
# plt.figure()
# Z = knn.predict(grid_e).reshape(xx.shape)
# plt.pcolormesh(grid_e[:, 0].reshape(xx.shape), grid_e[:, 1].reshape(xx.shape),
#             Z, cmap=cmap_light, alpha=.8)
# plt.scatter(X_embedded[:, 0], X_embedded[:, 1], c=y, cmap=cmap_bold, s=20, \
#                                                               edgecolor='k')
# plt.show()
