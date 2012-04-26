"""
===============================================================================
Multi-dimensional scaling - Reconstructing the map of the USA
===============================================================================

"""

# Author: Nelle Varoquaux <nelle.varoquaux@gmail.com>
# Licence: BSD

print __doc__
import numpy as np

from matplotlib import pyplot as plt
from matplotlib.collections import LineCollection

from sklearn import manifold
from sklearn.datasets import load_cities
from sklearn.metrics import euclidean_distances

cities_dataset = load_cities()
similarities = cities_dataset.data[:10, :10] * 1e-5

old_stress, best_pos = None, None
for n in range(10):
    mds = manifold.MDS(out_dim=2, n_init=6, max_iter=3000, n_jobs=2,
                    verbose=1, eps=1e-9)
    pos = mds.fit(similarities).positions_
    mds = manifold.MDS(out_dim=2, metric=False, n_init=1, max_iter=3000, n_jobs=2,
                    verbose=1, eps=1e-19)
    pos = mds.fit(similarities, init=pos).positions_
    stress = mds.stress_

    if old_stress is None or old_stress > stress:
        old_stress = stress
        best_pos = pos.copy()



fig = plt.figure(1)
distances = euclidean_distances(pos)
plt.scatter(distances, similarities, cmap=plt.cm.hot_r)

fig = plt.figure(2)
ax = plt.axes([0., 0., 1., 1.])
plt.scatter(pos[:, 0], pos[:, 1])

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

for index, (label, (x, y)) in enumerate(zip(cities_dataset.header, pos)):
    dx = x - pos[:, 0]
    dx[index] = 1
    dy = y - pos[:, 1]
    dy[index] = 1
    this_dx = dx[np.argmin(np.abs(dy))]
    this_dy = dy[np.argmin(np.abs(dx))]
    if this_dx > 0:
        horizontalalignment = 'left'
        x = x + .002
    else:
        horizontalalignment = 'right'
        x = x - .002
    if this_dy > 0:
        verticalalignment = 'bottom'
        y = y + .002
    else:
        verticalalignment = 'top'
        y = y - .002
    plt.text(x, y, label, size=10,
            horizontalalignment=horizontalalignment,
            verticalalignment=verticalalignment,
            bbox=dict(facecolor='w',
                      alpha=.6))

plt.title("Map of the US inferred from travel distances")
plt.show()
