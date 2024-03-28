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

no_of_iter = 20
no_of_points = 5000
no_of_clusters = 4
range_of_values = 100000
linkage_mats = [np.empty((0, no_of_clusters), int),
                np.empty((0, no_of_clusters), int),
                np.empty((0, no_of_clusters), int)]
linkage_type = ['ward', 'average', 'complete']
figure = [0, 0, 0]
final_plot = [0, 0, 0]
cluster_no = list(range(1, no_of_iter + 1))
np.random.seed(0)

for i in range(no_of_iter):
    X = np.random.randint(range_of_values, size=(no_of_points, 2))

    for j in range(3):
        clustering = AgglomerativeClustering(linkage=linkage_type[j],
                                             n_clusters=no_of_clusters)
        t0 = time()
        clustering.fit(X)
        print("%s : %.2fs" % (linkage_type[j], time() - t0))
        Y = np.zeros((no_of_clusters, ), dtype=np.int)

        for k in range(no_of_points):
            Y[clustering.labels_[k]] += 1

        Y = Y[np.argsort(Y)][::-1]
        linkage_mats[j] = np.append(linkage_mats[j],
                                    np.reshape(Y, (1, no_of_clusters)), axis=0)

for i in range(3):
    linkage_mats[i] = linkage_mats[i][np.argsort(linkage_mats[i][:, 0])]
    linkage_mats[i] = linkage_mats[i][:: -1, :]
    figure[i] = plt.figure()
    final_plot[i] = figure[i].add_subplot(111)
    final_plot[i].stackplot(cluster_no, linkage_mats[i].transpose())
    final_plot[i].set_title(linkage_type[i])
    final_plot[i].set_ylabel('Number of Points')
    final_plot[i].set_xlabel('Cluster Number')

plt.show()
