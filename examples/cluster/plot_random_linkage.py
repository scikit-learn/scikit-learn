"""
======================================================================================
Various Agglomerative Clustering on a 2D embedding of digits
======================================================================================
An illustration of various linkage option for agglomerative clustering on
a randomly generated dataset.

What this example shows us is the behavior "rich getting richer" of
agglomerative clustering that tends to create uneven cluster sizes.

Graphs help to visualize distribution of cluster sizes for each linkage type
and compare where behavior is more pronounced between average, ward & complete
linkage strategy. We can see that Ward creates more equally sized clusters.
Using random data helps to strengthen this notion.
"""

# Authors: Cheral Khandediya

from time import time
import numpy as np
from matplotlib import pyplot as plt
from sklearn.cluster import AgglomerativeClustering
import pandas as pd

ward_mat = []
average_mat = []
complete_mat = []


for j in range(50):

    array = np.random.randint(100000, size=(5000, 2))

    for linkage in ('ward', 'average', 'complete'):
        clustering = AgglomerativeClustering(linkage=linkage, n_clusters=4)
        t0 = time()
        clustering.fit(array)
        print("%s : %.2fs" % (linkage, time() - t0))

        arr = [0, 0, 0, 0]
        for i in range(5000):
            arr[clustering.labels_[i]] += 1
        arr.sort()
        mat = linkage + "_mat"
        eval(mat).append(arr)

ward_mat.sort()
complete_mat.sort()
average_mat.sort()

df_ward = pd.DataFrame(ward_mat, columns=list("1234"))
df_avg = pd.DataFrame(average_mat, columns=list("1234"))
df_comp = pd.DataFrame(complete_mat, columns=list("1234"))


# plots
df_comp.plot(kind='area')
df_avg.plot(kind='area')
df_ward.plot(kind='area')

plt.show()
