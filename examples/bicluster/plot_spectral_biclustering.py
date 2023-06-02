"""
=============================================
A demo of the Spectral Biclustering algorithm
=============================================

This example demonstrates how to generate a checkerboard dataset and bicluster
it using the Spectral Biclustering algorithm.

The data is generated, then shuffled and passed to the Spectral Biclustering
algorithm. The rows and columns of the shuffled matrix are then rearranged to
plot the biclusters found by the algorithm.
"""

# Author: Kemal Eren <kemal@kemaleren.com>
# License: BSD 3 clause

# %%
# Generate sample data
# --------------------
# We generate the sample data using the `make_checkerboard` function. Each pixel
# within `shape=(300, 300)` represents with it's color a value from a uniform
# distribution. The noise is added from a normal distribution, where the value
# chosen for `noise` is the standard deviation.
#
# As you can see, the data is distributed over 12 cluster cells and is
# relatively well distinguishable.

from sklearn.datasets import make_checkerboard
from matplotlib import pyplot as plt

n_clusters = (4, 3)
data, rows, columns = make_checkerboard(
    shape=(300, 300), n_clusters=n_clusters, noise=10, shuffle=False, random_state=1
)

plt.matshow(data, cmap=plt.cm.Blues)
plt.title("Original dataset")
plt.show()

# %%
# We're now going to shuffle the data. The goal is to reconstruct it afterwards
# using `SpectralBiclustering`.

import numpy as np

# Creating lists of shuffled row and column indices
rng = np.random.RandomState(0)
row_idx = rng.permutation(data.shape[0])
col_idx = rng.permutation(data.shape[1])

# %%
# We will redefine the shuffled data and plot it.
data = data[row_idx][:, col_idx]

plt.matshow(data, cmap=plt.cm.Blues)
plt.title("Shuffled dataset")
plt.show()

# %%
# Fitting SpectralBiclustering
# ----------------------------
# We fit the model and compare the obtained clusters with the ground truth
# biclusters (our artifical dataset). Note that when instanciating the
# model we specify the same number of clusters we used to create the dataset
# (`n_clusters = (4, 3)`), which will contribute to good results.

from sklearn.cluster import SpectralBiclustering
from sklearn.metrics import consensus_score

model = SpectralBiclustering(n_clusters=n_clusters, method="log", random_state=0)
model.fit(data)

# Calculating the similarity of two sets of biclusters
score = consensus_score(model.biclusters_, (rows[:, row_idx], columns[:, col_idx]))
print("consensus score: {:.1f}".format(score))

# %%
# The score is between 0 and 1. It shows the quality of the biclustering
# results.

# %%
# Plotting results
# ----------------
# Now, we rearrange the data based on the row and column labels assigned by the
# `SpectralBiclustering` model in ascending order and plot again. The
# `row_labels_` range from 0 to 3, while the `column_labels_` range from 0 to 2,
# representing a total of 4 clusters per row and 3 clusters per column.

fit_data = data[np.argsort(model.row_labels_)]
fit_data = fit_data[:, np.argsort(model.column_labels_)]

plt.matshow(fit_data, cmap=plt.cm.Blues)
plt.title("After biclustering; rearranged to show biclusters")
plt.show()


# %%
# As a last step, we want to demonstrate the relationships between the row
# and column labels assigned by the model. Therefore, we create a grid with
# `np.outer()`, which takes the sorted `row_labels_` and `column_labels_` and
# adds 1 to each to ensure that the labels start from 1 instead of 0 for better
# visualization.
#
# The outer product of the row and column label vectors show a representation
# of the checkerboard structure, where different combinations of row and column
# labels are represented by different shades of blue.

plt.matshow(
    np.outer(np.sort(model.row_labels_) + 1, np.sort(model.column_labels_) + 1),
    cmap=plt.cm.Blues,
)
plt.title("Checkerboard structure of rearranged data")
plt.show()
