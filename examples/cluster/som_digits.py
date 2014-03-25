# -*- coding: utf-8 -*-
"""
=======================================================================
A demo of Self-Organising Map and KMeans on the handwritten digits data
=======================================================================

Comparing various SOM and Kmeans clustering on the handwritten digits data
with the Caliński-Harabasz criterion (Caliński & Harabasz 1974).

"""

from __future__ import division

print(__doc__)

from time import time
import numpy as np

from sklearn.cluster import KMeans
from sklearn.cluster import SelfOrganizingMap
from sklearn.datasets import load_digits
from sklearn.preprocessing import scale

np.random.seed(42)


def Calinski_Harabasz_criterion(X, labels, centroids):
    '''Caliński-Harabasz criterion (pseudo F) statistic

    Ratio between the Between-group Sum of Squares (BGSS) and the
    Within-group Sum of Squares (WGSS), which acts similarly to an F statistic.
    The use is suggested by (Caliński & Harabasz 1974, p10) as an informal
    indicator of the best number of clusters (k) for a given sample (size n).

    C-H criterion = [BGSS/(k - 1)] / [WGSS/(n - k)]

    - Caliński, T. & Harabasz, J., 1974. A dendrite method for cluster
      analysis. Communications in Statistics, 3(1), pp.1-27.
      doi://10.1080/03610927408827101
    '''
    mean = np.mean(X, axis=0)
    BGSS = np.sum([(k - mean)**2 for k in centroids])
    WGSS = np.sum([(x - centroids[labels[i]])**2 for i, x in enumerate(X)])
    k = len(centroids)
    n = len(X)
    return (BGSS / (k-1)) / (WGSS / (n-k))


# Load dataset

digits = load_digits()
data = scale(digits.data)
n_samples, n_features = data.shape
n_digits = len(np.unique(digits.target))

print("Digits dataset")
print("n_digits   : %d" % n_digits)
print("n_features : %d" % n_features)
print("n_samples  : %d" % n_samples)
print()

# Digits dataset clustering using Self-Organizing Map

print("Self-Organizing Map ")
t0 = time()
grid_width = 4
som = SelfOrganizingMap(affinity=(grid_width, grid_width),
                        n_iterations=n_samples*5, learning_rate=1)
som.fit(data)
print("done in %0.3fs" % (time() - t0))
print()

F = Calinski_Harabasz_criterion(data, som.labels_, som.centers_)
print('Caliński-Harabasz criterion (pseudo F): %0.2f | %0.2f%%' %
      (F, 100 * (F / (1 + F))))
print()

# Digits dataset clustering using Kmeans

print("KMeans ")
t0 = time()
km = KMeans(init='k-means++', n_clusters=grid_width**2, n_init=10)
km.fit(data)
print("done in %0.3fs" % (time() - t0))
print()

F = Calinski_Harabasz_criterion(data, km.labels_, km.cluster_centers_)
print('Caliński-Harabasz criterion (pseudo F): %0.2f | %0.2f%%' %
      (F, 100 * (F / (1 + F))))
