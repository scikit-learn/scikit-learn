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


digits = load_digits()
X, y = digits.data, digits.target
tsne = TSNE(n_components=2, perplexity=40, learning_rate=500, verbose=2,
            random_state=0)
X_embedded = tsne.fit_transform(X)

plt.figure(figsize=(6, 5))
plt.title("Digits dataset, trustworthiness = %.3f"
          % tsne.score(X, n_neighbors=12))
plt.scatter(X_embedded[:, 0], X_embedded[:, 1], c=y)
plt.xticks(())
plt.yticks(())
plt.show()
