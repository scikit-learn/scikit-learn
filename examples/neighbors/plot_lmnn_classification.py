"""
============================================
Large Margin Nearest Neighbor Classification
============================================

Sample usage of Large Margin Nearest Neighbor classification.
It will plot the decision boundaries for each class.
"""

# Author: John Chiotellis <johnyc.code@gmail.com>
# License: BSD 3 clause

import numpy as np
import matplotlib.pyplot as plt
from matplotlib.colors import ListedColormap
from sklearn import datasets
from sklearn.model_selection import train_test_split
from sklearn.neighbors import LargeMarginNearestNeighbor, KNeighborsClassifier


print(__doc__)

n_neighbors = 3

# import some data to play with
iris = datasets.load_iris()

# we only take the first two features. We could avoid this ugly
# slicing by using a two-dim dataset
X = iris.data[:, :2]
y = iris.target

X_train, X_test, y_train, y_test = \
    train_test_split(X, y, stratify=y, test_size=0.7, random_state=42)

h = .02  # step size in the mesh

# Create color maps
cmap_light = ListedColormap(['#FFAAAA', '#AAFFAA', '#AAAAFF'])
cmap_bold = ListedColormap(['#FF0000', '#00FF00', '#0000FF'])


def plot_feature_space(ax, lmnn=None):
    # we create an instance of Neighbours Classifier and fit the data.
    clf = KNeighborsClassifier(n_neighbors)

    # Plot the decision boundary. For that, we will assign a color to each
    # point in the mesh [x_min, x_max]x[y_min, y_max].
    x_min, x_max = X[:, 0].min() - 1, X[:, 0].max() + 1
    y_min, y_max = X[:, 1].min() - 1, X[:, 1].max() + 1
    xx, yy = np.meshgrid(np.arange(x_min, x_max, h),
                         np.arange(y_min, y_max, h))

    if lmnn is None:
        clf.fit(X_train, y_train)
        space = np.c_[xx.ravel(), yy.ravel()]
        score = clf.score(X_test, y_test)
        method = 'K-Nearest Neighbors'
    else:
        clf.fit(lmnn.transform(X_train), y_train)
        space = lmnn.transform(np.c_[xx.ravel(), yy.ravel()])
        score = clf.score(lmnn.transform(X_test), y_test)
        method = 'Large Margin Nearest Neighbor'

    Z = clf.predict(space)

    # Put the result into a color plot
    Z = Z.reshape(xx.shape)
    ax.contourf(xx, yy, Z, cmap=cmap_light, alpha=.8)

    # Plot also the training points
    ax.scatter(X[:, 0], X[:, 1], c=y, cmap=cmap_bold, edgecolor='k', s=20)
    ax.set_xlim(xx.min(), xx.max())
    ax.set_ylim(yy.min(), yy.max())
    ax.set_title("{} (k = {})".format(method, n_neighbors))
    ax.text(xx.max() - .3, yy.min() + .3, ('%.2f' % score).lstrip('0'),
            size=15, horizontalalignment='right')


fig = plt.figure(figsize=(15, 6))

ax = fig.add_subplot(121)
plot_feature_space(ax)

lmnn = LargeMarginNearestNeighbor(n_neighbors=n_neighbors, random_state=42)
lmnn.fit(X_train, y_train)

ax = fig.add_subplot(122)
plot_feature_space(ax, lmnn=lmnn)

plt.show()
