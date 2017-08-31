"""
======================================================
Large Margin Nearest Neighbor Dimensionality Reduction
======================================================

Sample usage of Large Margin Nearest Neighbor for dimensionality reduction.
"""

# Author: John Chiotellis <johnyc.code@gmail.com>
#
# License: BSD 3 clause

import numpy as np
import matplotlib.pyplot as plt
from sklearn import neighbors, datasets, decomposition, model_selection
from sklearn.discriminant_analysis import LinearDiscriminantAnalysis


print(__doc__)

n_neighbors = 3
random_state = 0

faces = datasets.fetch_olivetti_faces()
X, y = faces.data, faces.target
X_train, X_test, y_train, y_test = model_selection.\
    train_test_split(X, y, test_size=0.5, stratify=y, random_state=random_state)

dim = len(X[0])
n_classes = len(np.unique(y))
cmap_bold = None

# Put the result into a color plot
fig = plt.figure()

# Reduce dimensions with PCA
pca = decomposition.PCA(n_components=2, random_state=random_state)
X_pca = pca.fit_transform(np.concatenate((X_train, X_test)))

# Compute nearest neighbor accuracy
knn = neighbors.KNeighborsClassifier(n_neighbors=n_neighbors)
knn.fit(X_pca[:len(X_train)], y_train)
acc_pca = knn.score(X_pca[len(X_train):], y_test)

# Plot the points after PCA
fig.add_subplot(131)
plt.scatter(X_pca[:, 0], X_pca[:, 1], c=y, cmap=cmap_bold)
plt.title("PCA (3-nn test acc. = {:5.2f}%)".format(acc_pca * 100))


# Reduce dimensions with LinearDiscriminantAnalysis
lda = LinearDiscriminantAnalysis(n_components=2)
lda = lda.fit(X_train, y_train)
LX = lda.transform(X)

# Compute nearest neighbor accuracy
knn.fit(lda.transform(X_train), y_train)
acc_lda = knn.score(lda.transform(X_test), y_test)

# Plot the points after transformation with LinearDiscriminantAnalysis
fig.add_subplot(132)
plt.scatter(LX[:, 0], LX[:, 1], c=y, cmap=cmap_bold)
plt.title("LDA (test acc. = {:5.2f}%)".format(acc_lda * 100))

# Reduce dimensions with LargeMarginNearestNeighbor
lmnn = neighbors.LargeMarginNearestNeighbor(n_neighbors=n_neighbors,
                                            n_features_out=2,
                                            verbose=1,
                                            max_iter=30,
                                            random_state=random_state)
lmnn = lmnn.fit(X_train, y_train)
LX = lmnn.transform(X)

# Compute nearest neighbor accuracy
knn.fit(lmnn.transform(X_train), y_train)
acc_lmnn = knn.score(lmnn.transform(X_test), y_test)

# Plot the points after transformation with LargeMarginNearestNeighbor
fig.add_subplot(133)
plt.scatter(LX[:, 0], LX[:, 1], c=y, cmap=cmap_bold)
plt.title("LMNN (k = {}, test acc. = {:5.2f}%)".format(n_neighbors,
                                                       acc_lmnn * 100))

plt.suptitle('Dimensionality Reduction ({} dimensions, {} classes)'.format(
    dim, n_classes), fontweight='bold', fontsize=18)

plt.show()
