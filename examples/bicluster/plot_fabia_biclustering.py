"""
==========================================
A demo of the Fabia Biclustering algorithm
==========================================

This example demonstrates how to generate a dataset and bicluster it
using the the Fabia Biclustering algorithm.

The dataset is generated using the ``make_fabia_biclusters`` function, which
implants small and large biclusters in a matrix as well as a lot of background
noise. It then runs the Fabia Biclustering algorithm on this dataset,
and overlays the found biclusters on top of the original data.
"""
print(__doc__)

# Author: Thomas Unterthiner
# License: BSD 3 clause

import matplotlib.pyplot as plt
from sklearn.cluster import FabiaBiclustering
from sklearn.metrics import consensus_score
from sklearn.datasets.samples_generator import make_fabia_biclusters

# Create data
k, n, m = 8, 500, 600
(Xn, X, zc, lc) = make_fabia_biclusters((n, m), k, (5, 100), (10, 250), 2.0,
                                        0.2, 2.0, 1.0, 0.2, 2.0, 1.0, as_blocks=True,
                                        random_state=42)

fig, ax = plt.subplots(figsize=(10, 12))
im = plt.imshow(Xn, cmap=plt.cm.RdBu)
ax.set_title("Input data")
fig.colorbar(im, ax=ax)
fig.show()

model = FabiaBiclustering(n_clusters=9, alpha=0.05, scale=False, random_state=42, thresZ=2.0)
model.fit(Xn)

# Show data with found biclusters overlayed
fig, ax = plt.subplots(figsize=(10, 12))
r = np.dot(model.Z_, model.L_)
im = ax.imshow(r, cmap=plt.cm.RdBu)
overlay = np.zeros_like(r)
for bic in range(model.n_clusters):
    for i in np.where(model.biclusters_[0][bic])[0]:
        overlay[i, model.biclusters_[1][bic]] = 1
ax.imshow(overlay, alpha=0.5, cmap=plt.cm.Greens)
ax.set_title("Overlayed biclusters")
fig.show()
