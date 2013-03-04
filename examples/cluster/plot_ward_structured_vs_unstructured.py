"""
===========================================================
Hierarchical clustering: structured vs unstructured ward
===========================================================

Example builds a swiss roll dataset and runs
:ref:`hierarchical_clustering` on their position.

In a first step, the hierarchical clustering without connectivity
constraints on structure, solely based on distance, whereas in a second
step clustering restricted to the k-Nearest Neighbors graph: it's a
hierarchical clustering with structure prior.

Some of the clusters learned without connectivity constraints do not
respect the structure of the swiss roll and extend across different folds of
the manifolds. On the opposite, when opposing connectivity constraints,
the clusters form a nice parcellation of the swiss roll.
"""

# Authors : Vincent Michel, 2010
#           Alexandre Gramfort, 2010
#           Gael Varoquaux, 2010
# License: BSD

print(__doc__)

import time as time
import numpy as np
import pylab as pl
import mpl_toolkits.mplot3d.axes3d as p3
from sklearn.cluster import Ward
from sklearn.datasets.samples_generator import make_swiss_roll

###############################################################################
# Generate data (swiss roll dataset)
n_samples = 1000
noise = 0.05
X, _ = make_swiss_roll(n_samples, noise)
# Make it thinner
X[:, 1] *= .5

###############################################################################
# Compute clustering
print("Compute unstructured hierarchical clustering...")
st = time.time()
ward = Ward(n_clusters=6).fit(X)
class_ = ward.classes_
print("Elapsed time: ", time.time() - st)
print("Number of points: ", class_.size)

###############################################################################
# Plot result
fig = pl.figure()
ax = p3.Axes3D(fig)
ax.view_init(7, -80)
for l in np.unique(class_):
    ax.plot3D(X[class_ == l, 0], X[class_ == l, 1], X[class_ == l, 2],
              'o', color=pl.cm.jet(np.float(l) / np.max(class_ + 1)))
pl.title('Without connectivity constraints')


###############################################################################
# Define the structure A of the data. Here a 10 nearest neighbors
from sklearn.neighbors import kneighbors_graph
connectivity = kneighbors_graph(X, n_neighbors=10)

###############################################################################
# Compute clustering
print("Compute structured hierarchical clustering...")
st = time.time()
ward = Ward(n_clusters=6, connectivity=connectivity).fit(X)
class_ = ward.classes_
print("Elapsed time: ", time.time() - st)
print("Number of points: ", class_.size)

###############################################################################
# Plot result
fig = pl.figure()
ax = p3.Axes3D(fig)
ax.view_init(7, -80)
for l in np.unique(class_):
    ax.plot3D(X[class_ == l, 0], X[class_ == l, 1], X[class_ == l, 2],
              'o', color=pl.cm.jet(float(l) / np.max(class_ + 1)))
pl.title('With connectivity constraints')

pl.show()
