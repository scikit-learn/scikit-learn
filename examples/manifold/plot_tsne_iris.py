"""
===================================
t-SNE Visualization of Iris Dataset
===================================

The iris dataset will be split randomly in a training set and a test set. We
use the training set to learn a two-dimensional embedding of the training data
which is displayed on the left side. It is possible to generalize and map the
test data to the embedded space as well. This is shown on the right side of
the plot.
"""

import numpy as np
import matplotlib.pyplot as plt
from sklearn.manifold import TSNE
from sklearn.datasets import load_iris
from sklearn.cross_validation import train_test_split


iris = load_iris()
X, y = iris.data, iris.target
X_train, X_test, y_train, y_test = train_test_split(X, y, test_size=0.5,
                                                    random_state=0)
tsne = TSNE(n_components=2, perplexity=20, learning_rate=100.0, random_state=0)
tsne.fit(X_train)
X_embedded_train = tsne.transform(X_train)
X_embedded_test = tsne.transform(X_test)

plt.figure(figsize=(12, 5))
plt.suptitle("Two-dimensional embedding of Iris dataset with t-SNE")
ax = plt.subplot(121)
ax.set_title("Training data, trustworthiness = %.3f" % tsne.score(X_train))
p1 = ax.scatter(X_embedded_train[y_train == 0, 0],
                X_embedded_train[y_train == 0, 1], c="r")
p2 = ax.scatter(X_embedded_train[y_train == 1, 0],
                X_embedded_train[y_train == 1, 1], c="g")
p3 = ax.scatter(X_embedded_train[y_train == 2, 0],
                X_embedded_train[y_train == 2, 1], c="b")
ax.legend((p1, p2, p3), iris.target_names, loc="center left")
plt.setp(ax, xticks=(), yticks=())
ax = plt.subplot(122)
ax.set_title("Test data, trustworthiness = %.3f" % tsne.score(X_test))
p1 = ax.scatter(X_embedded_test[y_test == 0, 0],
                X_embedded_test[y_test == 0, 1], c="r")
p2 = ax.scatter(X_embedded_test[y_test == 1, 0],
                X_embedded_test[y_test == 1, 1], c="g")
p3 = ax.scatter(X_embedded_test[y_test == 2, 0],
                X_embedded_test[y_test == 2, 1], c="b")
ax.legend((p1, p2, p3), iris.target_names, loc="center left")
plt.setp(ax, xticks=(), yticks=())
plt.show()
