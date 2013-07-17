"""
=====================================================
Applying Spectral Biclustering to a generated dataset
=====================================================

This example demonstrates how to generate a checkerboard dataset for
the Spectral Biclustering algorithm and fit it.

The data is generated with the `make_checkerboard` function, then
shuffled and passed to the Spectral Biclustering algorithm. The rows
and columns of the shuffled matrix are rearranged to show the
biclusters found by the algorithm.

The outer product of the row and column label vectors shows a
representation of the checkerboard structure.

"""
print(__doc__)

# Author: Kemal Eren <kemal@kemaleren.com>
# License: BSD 3 clause

import numpy as np
import pylab as pl

from sklearn.datasets import make_checkerboard
from sklearn.cluster.bicluster import SpectralBiclustering

n_clusters = (4, 3)
data, row_labels, column_labels = make_checkerboard(
    shape=(300, 300), n_clusters=n_clusters, noise=10,
    shuffle=True, random_state=0)
model = SpectralBiclustering(n_clusters=n_clusters, method='log',
                             random_state=0)
model.fit(data)

fit_data = data[np.argsort(model.row_labels_)]
fit_data = fit_data[:, np.argsort(model.column_labels_)]

pl.matshow(data[np.argsort(row_labels)]
               [:, np.argsort(column_labels)],
           cmap=pl.cm.Blues)
pl.title("Original dataset")

pl.matshow(data, cmap=pl.cm.Blues)
pl.title("Shuffled dataset")

pl.matshow(fit_data, cmap=pl.cm.Blues)
pl.title("After biclustering; rearranged to show biclusters")

pl.matshow(np.outer(np.sort(model.row_labels_) + 1,
                    np.sort(model.column_labels_) + 1),
           cmap=pl.cm.Blues)
pl.title("Checkerboard structure of rearranged data")

pl.show()
