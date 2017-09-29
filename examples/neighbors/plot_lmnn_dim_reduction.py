"""
======================================================
Large Margin Nearest Neighbor Dimensionality Reduction
======================================================

Sample usage of Large Margin Nearest Neighbor for dimensionality reduction.
"""

# Author: John Chiotellis <johnyc.code@gmail.com>
# License: BSD 3 clause

import numpy as np
import matplotlib.pyplot as plt
from sklearn import datasets
from sklearn.model_selection import train_test_split
from sklearn.decomposition import PCA
from sklearn.discriminant_analysis import LinearDiscriminantAnalysis
from sklearn.neighbors import LargeMarginNearestNeighbor, KNeighborsClassifier


print(__doc__)

n_neighbors = 3
random_state = 0

# Load Olivetti Faces dataset
faces = datasets.fetch_olivetti_faces()
X, y = faces.data, faces.target

# Split into train/test
X_train, X_test, y_train, y_test = \
    train_test_split(X, y, test_size=0.5, stratify=y,
                     random_state=random_state)

dim = len(X[0])
n_classes = len(np.unique(y))

# Reduce dimension to 2 with PCA
pca = PCA(n_components=2, random_state=random_state)

# Reduce dimension to 2 with LinearDiscriminantAnalysis
lda = LinearDiscriminantAnalysis(n_components=2)

# Reduce dimension to 2 with LargeMarginNearestNeighbor
lmnn = LargeMarginNearestNeighbor(n_neighbors=n_neighbors, n_features_out=2,
                                  max_iter=20, verbose=1,
                                  random_state=random_state)

# Use a nearest neighbor classifier to evaluate the methods
knn = KNeighborsClassifier(n_neighbors=n_neighbors)

# Make a list of the methods to be compared
dim_reduction_methods = [('PCA', pca), ('LDA', lda), ('LMNN', lmnn)]

for i, (name, method) in enumerate(dim_reduction_methods):
    plt.figure()

    # Fit the method's model
    method.fit(X_train, y_train)

    # Fit a nearest neighbor classifier on the embedded training set
    knn.fit(method.transform(X_train), y_train)

    # Compute the nearest neighbor accuracy on the embedded test set
    acc_knn = knn.score(method.transform(X_test), y_test)

    # Embed the data set in 2 dimensions using the method
    X_embedded = method.transform(X)

    # Plot the embedding and show the evaluation score
    plt.scatter(X_embedded[:, 0], X_embedded[:, 1], c=y)
    plt.title("{}, KNN (k={})".format(name, n_neighbors))
    plt.text(0.9, 0.1, '{:.2f}'.format(acc_knn), size=15,
             ha='center', va='center', transform=plt.gca().transAxes)

plt.show()
