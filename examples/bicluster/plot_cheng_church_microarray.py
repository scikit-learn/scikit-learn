"""
==================================================
Biclustering microarray data with Cheng and Church
==================================================

This example is a replication of an experiment from the original Cheng
and Church paper. It demonstrates how to use this model to bicluster a
gene expression microarray dataset. Each row represents a different
gene, and each column represents a tissue sample from a patient with
lymphoma. The larger the value of ``data[i, j]``, the more active gene
``i`` in sample ``j``. Biclustering this data with Cheng and Church
finds subsets of samples with similar expression profiles in a subset
of genes. The goal of this kind of analysis is often to find sets of
genes that may be somehow related. For instance, lymphoma may cause
some genes that are otherwise unrelated to become highly expressed or
supressed.

The gene microarray data is downloaded from the paper's supplementary
information webpage, parsed into a NumPy array, and clustered with
Cheng and Church. The bicluster is then visualized by a parallel
coordinate plot of its rows. Biclustering is performed with almost the
same parameters as in the original experiment, except the mean squared
residue threshold is lowered to make the bicluster visually simpler.

"""
from __future__ import print_function

print(__doc__)

from time import time
import urllib

import numpy as np
from matplotlib import pyplot as plt

from sklearn.cluster.bicluster import ChengChurch

# get data
url = "http://arep.med.harvard.edu/biclustering/lymphoma.matrix"
lines = urllib.urlopen(url).read().strip().split('\n')
# insert a space before all negative signs
lines = list(' -'.join(line.split('-')).split(' ') for line in lines)
lines = list(list(int(i) for i in line if i) for line in lines)
data = np.array(lines)

# replace missing values, just as in the paper
generator = np.random.RandomState(0)
idx = np.where(data == 999)
data[idx] = generator.randint(-800, 801, len(idx[0]))

# cluster with similar parameters as original paper
model = ChengChurch(n_clusters=1, max_msr=200,
                    deletion_threshold=1.2, inverse_rows=True,
                    random_state=0)
print("Biclustering...")
start_time = time()
model.fit(data)
print("Done in {:.2f}s.".format(time() - start_time))
bicluster = model.get_submatrix(0, data)

msr = lambda a: (np.power(a - a.mean(axis=1, keepdims=True) -
                          a.mean(axis=0) + a.mean(), 2).mean())

print("Bicluster MSR: {}".format(msr(bicluster)))

n_cols = bicluster.shape[1]
for row in bicluster:
    plt.plot(np.arange(n_cols), row)
plt.title('Parallel coordinates of the bicluster rows')
plt.xlabel('column number')
plt.ylabel('expression level')
plt.show()
