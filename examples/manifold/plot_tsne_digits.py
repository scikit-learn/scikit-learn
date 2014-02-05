"""
=====================================
t-SNE Visualization of Digits Dataset
=====================================

TODO
"""

import numpy as np
import matplotlib.pyplot as plt
from sklearn.manifold import TSNE
from sklearn.datasets import load_digits


if __name__ == "__main__":
    digits = load_digits()
    X, y = digits.data, digits.target
    tsne = TSNE(n_components=2, optimizer="gradient_descent", perplexity=40,
                n_iter=1000, verbose=1, random_state=0)
    X_embedded = tsne.fit_transform(X)

    plt.figure(figsize=(6, 5))
    plt.title("Digits dataset, trustworthiness = %.3f"
              % tsne.score(X, n_neighbors=12))
    plt.scatter(X_embedded[:, 0], X_embedded[:, 1], c=y)
    plt.xticks(())
    plt.yticks(())
    plt.show()
