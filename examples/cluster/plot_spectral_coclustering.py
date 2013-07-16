"""
======================================================
Applying Spectral Co-Clustering to a generated dataset
======================================================

This example demonstrates how to generate a dataset for the Spectral
Co-Clustering algorithm and fit it.

The dataset is generated using the `make_biclusters` function, which
creates a matrix of small values and implants bicluster with large
values. The rows and columns are then shuffled and passed to the
Spectral Co-Clustering algorithm. Rearranging the shuffled matrix to
make biclusters contiguous shows how accurately the algorithm found
the biclusters.

"""
print(__doc__)

# Author: Kemal Eren <kemal@kemaleren.com>
# License: BSD 3 clause

import numpy as np
import pylab as pl

from sklearn.datasets import make_biclusters
from sklearn.cluster.bicluster import SpectralCoclustering

data, row_labels, column_labels = make_biclusters(
    shape=(300, 300), n_clusters=5, noise=10,
    shuffle=True, random_state=0)

model = SpectralCoclustering(n_clusters=5)
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

pl.show()
