"""
========================================
A demo of the Cheng and Church algorithm
========================================

This example demonstrates how to generate a dataset and bicluster it
using Cheng and Church.

The data is generated with the ``make_msr_biclusters`` function, then
shuffled and passed to Cheng and Church. The rows and columns of the
shuffled matrix are rearranged to show the biclusters found by the
algorithm.

"""

print(__doc__)

# Author: Kemal Eren <kemal@kemaleren.com>
# License: BSD 3 clause

from matplotlib import pyplot as plt

from sklearn.datasets import make_msr_biclusters
from sklearn.datasets import samples_generator as sg
from sklearn.cluster.bicluster import ChengChurch
from sklearn.metrics import consensus_score

data, rows, columns = make_msr_biclusters(shape=(100, 100),
                                          n_clusters=3, noise=10,
                                          shuffle=False,
                                          random_state=0)

plt.matshow(data, cmap=plt.cm.Blues)
plt.title("Original dataset")

data, row_idx, col_idx = sg._shuffle(data, random_state=0)

plt.matshow(data, cmap=plt.cm.Blues)
plt.title("Shuffled dataset")

plt.figure()
n_cols = data.shape[1]
for row in data:
    plt.plot(range(n_cols), row)
plt.title('Parallel coordinates of shuffled data')
plt.xlabel('column numbers')
plt.ylabel('value')

model = ChengChurch(n_clusters=3, max_msr=100, random_state=0)
model.fit(data)
score = consensus_score(model.biclusters_,
                        (rows[:, row_idx], columns[:, col_idx]))

print "consensus score: {:.1f}".format(score)


plt.figure()
bicluster = model.get_submatrix(0, data)
n_cols = bicluster.shape[1]
for row in bicluster:
    plt.plot(range(n_cols), row)
plt.title('Parallel coordinates of first bicluster')
plt.xlabel('column numbers')
plt.ylabel('value')
plt.xlim(0, n_cols)

plt.show()
