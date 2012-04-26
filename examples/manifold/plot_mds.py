"""
=========================
Multi-dimensional scaling
=========================

"""

# Author: Nelle Varoquaux <nelle.varoquaux@gmail.com>
# Licence: BSD

print __doc__
import numpy as np

from matplotlib import pyplot as plt
from matplotlib.collections import LineCollection

from sklearn import manifold
from sklearn.metrics import euclidean_distances
from sklearn.decomposition import PCA

n_samples = 20
X_true = np.random.randint(0, 20, 2 * n_samples)
X_true = X_true.reshape((n_samples, 2))
# Center the data
X_true -= X_true.mean()

similarities = euclidean_distances(X_true)

# Add noise to the similarities
noise = np.random.rand(n_samples, n_samples)
noise += noise.T
noise[np.arange(noise.shape[0]), np.arange(noise.shape[0])] = 0
similarities += noise

mds = manifold.MDS(out_dim=2, max_iter=3000, n_jobs=2,
                   eps=1e-9)
pos = mds.fit(similarities).positions_

nmds = manifold.MDS(out_dim=2, metric=False,
                    max_iter=3000, n_jobs=2,
                    eps=1e-9)
npos = mds.fit(similarities).positions_

# Rotate the data
clf = PCA(n_components=3)
X_true = clf.fit_transform(X_true)

pos = clf.fit_transform(pos)

npos = clf.fit_transform(pos)

fig = plt.figure(1)
ax = plt.axes([0., 0., 1., 1.])

plt.scatter(X_true[:, 0] + 0.2, X_true[:, 1] + 0.2, c='r', s=10)
plt.scatter(pos[:, 0] + 0.2, pos[:, 1] + 0.2, s=10, c='g')
plt.scatter(pos[:, 0] - 0.2, pos[:, 1] - 0.2, s=10, c='b')
plt.legend(('True position', 'MDS', 'NMDS'))

similarities = similarities.max() / similarities * 100
similarities[np.isinf(similarities)] = 0

# Plot the edges
start_idx, end_idx = np.where(pos)
#a sequence of (*line0*, *line1*, *line2*), where::
#            linen = (x0, y0), (x1, y1), ... (xm, ym)
segments = [[pos[i, :], pos[j, :]]
            for i in range(len(pos)) for j in range(len(pos))]
values = np.abs(similarities)
lc = LineCollection(segments,
                    zorder=0, cmap=plt.cm.hot_r,
                    norm=plt.Normalize(0, values.max()))
lc.set_array(similarities.flatten())
lc.set_linewidths(0.5 * np.ones(len(segments)))
ax.add_collection(lc)

plt.show()
