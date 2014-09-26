# -*- coding: utf-8 -*-
"""
Created on Wed Sep 24 10:18:50 2014

@author: clemens
"""

from sklearn.datasets import make_blobs
from sklearn.lda import LDA
import matplotlib.pyplot as plt

import numpy as np


def plot_data(X, y):
    for label, marker, color in zip(np.unique(y), ("x", "o"), ("blue", "red")):
        plt.scatter(x=X[y == label, 0], y=X[y == label, 1], marker=marker, color=color)


def plot_boundary(clf):
    nx, ny = 200, 100
    x_min, x_max = plt.xlim()
    y_min, y_max = plt.ylim()
    xx, yy = np.meshgrid(np.linspace(x_min, x_max, nx),
                         np.linspace(y_min, y_max, ny))
    Z = clf.predict_proba(np.c_[xx.ravel(), yy.ravel()])
    Z = Z[:, 1].reshape(xx.shape)
    plt.contour(xx, yy, Z, [0.5], linewidths=2., colors='k')

# Generate data
X, y = make_blobs(n_samples=100, n_features=2, centers=[[-2, 5], [0, 3.5]])

clf1 = LDA(solver="svd", store_covariance=True)
clf2 = LDA(solver="lsqr", alpha=None)
clf3 = LDA(solver="eigen", alpha=None)
clf1.fit(X, y)
clf2.fit(X, y)
clf3.fit(X, y)

plt.figure(1)
plt.subplot(1, 3, 1)
plot_data(X, y)
plot_boundary(clf1)
plt.title("svd")

plt.subplot(1, 3, 2)
plot_data(X, y)
plot_boundary(clf2)
plt.title("lsqr")

plt.subplot(1, 3, 3)
plot_data(X, y)
plot_boundary(clf3)
plt.title("eigen")
plt.show()
