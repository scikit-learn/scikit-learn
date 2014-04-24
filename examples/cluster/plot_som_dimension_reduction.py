"""
=============================================================
A demo of SelfOrganisingMap used for dimensionality reduction
=============================================================

Example of SOM clustering used for non-linear dimensionality reduction. We
start with a PCA-based initialisation, and use a Self-Organising Map to fit
a non-linear pattern in the data. 

"""

print(__doc__)

import matplotlib.pyplot as plt
import numpy as np
from sklearn.cluster import SelfOrganizingMap
from sklearn.decomposition import PCA


def noisy_s(t):
    stdev = 0.2
    xy = (t + stdev * np.random.randn(),
          np.sin(t) + stdev * np.random.randn())
    return (xy)


def plot_line_2_points(x, y, xlims):
    """from http://stackoverflow.com/a/9153411/210945"""
    xrange = np.arange(xlims[0], xlims[1], 0.1)
    A = np.vstack([x, np.ones(len(x))]).T
    k, b = np.linalg.lstsq(A, y)[0]
    plt.plot(xrange, k*xrange + b, 'b')

# Data size
n_samples = 2000
data = np.array([noisy_s(np.pi * (2*t / n_samples - 1)) for t in range(n_samples)])

# SOM size
n_centers = 10

som = SelfOrganizingMap(adjacency=(n_centers, 1), n_iterations=1000,
                        init='pca', learning_rate=1)
som.fit(data)

pca = PCA()
pca.fit(data)
basis = np.array([[[-1,0], [1,0]],[[0,-1], [0,1]]])
pca_basis = pca.inverse_transform(basis)

plt.scatter(data[:, 0], data[:, 1], c=som.labels_, cmap=plt.cm.prism, alpha=0.3)
plt.plot(som.cluster_centers_[:, 0], som.cluster_centers_[:, 1], '-r', linewidth=2)
plt.plot(som.cluster_centers_[:, 0], som.cluster_centers_[:, 1], 'or')
for b in pca_basis:
    plot_line_2_points(b[0], b[1], [2,2])
plt.title('Organized Map')
F = plt.gcf()
F.set_size_inches((40, 20))
plt.show()
