"""
==================================================
Biclustering microarray data with Cheng and Church
==================================================

This example is a replication of an experiment from the original Cheng
and Church paper. The gene microarray data is downloaded from the
paper's supplementary information webpage and parsed into a NumPy
array. Cheng and Church is used with the same parameters as in the
original experiment to find 100 biclusters.

Output::

    Biclustering...
    Done in 36.28s.
    MSR of best bicluster: 853.14

"""
from __future__ import print_function

print(__doc__)

from time import time
import urllib

import numpy as np

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

# cluster with same parameters as original paper
model = ChengChurch(n_clusters=100, max_msr=1200,
                    deletion_threshold=1.2, inverse_rows=True,
                    random_state=0)
print("Biclustering...")
start_time = time()
model.fit(data)
print("Done in {:.2f}s.".format(time() - start_time))

# find smallest msr
msr = lambda a: (np.power(a - a.mean(axis=1, keepdims=True) -
                          a.mean(axis=0) + a.mean(), 2).mean())
min_msr = min(msr(model.get_submatrix(i, data)) for i in range(100))
print ("MSR of best bicluster: {:.2f}".format(min_msr))
