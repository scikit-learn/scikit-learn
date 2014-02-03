"""
===================================
t-SNE Visualization of Iris Dataset
===================================

TODO
"""

import numpy as np
import matplotlib.pyplot as plt
from sklearn.manifold import TSNE
from sklearn.datasets import load_iris


if __name__ == "__main__":
    iris = load_iris()
    X, y = iris.data, iris.target
    tsne = TSNE(n_components=2, optimizer="gradient_descent",
                perplexity=20, n_iter=1000, verbose=1,
                random_state=0)
    X_embedded = tsne.fit_transform(X)

    plt.title("Two-dimensional embedding of Iris dataset with t-SNE")
    ax = plt.subplot(111)
    p1 = ax.scatter(X_embedded[y == 0, 0], X_embedded[y == 0, 1], c="r")
    p2 = ax.scatter(X_embedded[y == 1, 0], X_embedded[y == 1, 1], c="g")
    p3 = ax.scatter(X_embedded[y == 2, 0], X_embedded[y == 2, 1], c="b")
    ax.legend((p1, p2, p3), iris.target_names)
    plt.setp(ax, xticks=(), yticks=())
    plt.show()
