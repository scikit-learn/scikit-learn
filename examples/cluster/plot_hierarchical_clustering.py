import time
import pylab as pl

from sklearn.cluster.hierarchical import HierarchicalLinkage
from sklearn.neighbors import radius_neighbors_graph
from sklearn import datasets

##############################################################################
# Generate sample data

N = 2000
# XXX: random_state=0 triggers a bug...
#X = datasets.make_moons(n_samples=N, noise=.05, random_state=0)[0]
X = datasets.make_swiss_roll(n_samples=N, noise=.5, random_state=42)[0]
X = X[:, [0, 2]]

knn_graph = radius_neighbors_graph(X, 1.5)

for n_clusters in (60, 20, 3):
    pl.figure(figsize=(12, 7))
    for connectivity in (None, knn_graph):
        for index, linkage in enumerate(('average', 'complete', 'ward')):
            pl.subplot(2, 3, (connectivity is None) * 3 + index + 1)
            model = HierarchicalLinkage(linkage=linkage,
                                        connectivity=connectivity,
                                        n_clusters=n_clusters)
            t0 = time.time()
            model.fit(X)
            elapsed_time = time.time() - t0
            pl.scatter(X[:, 0], X[:, 1], c=model.labels_, cmap=pl.cm.spectral)
            pl.title('linkage=%s, connectivity=%r \n(time %.2fs)' % (linkage,
                     connectivity is not None, elapsed_time),
                     fontdict=dict(verticalalignment='top'))
            pl.axis('off')

        pl.subplots_adjust(bottom=0, top=.9, hspace=0, wspace=0,
                           left=0, right=1)
        pl.suptitle('n_cluster=%i' % n_clusters, size=15)


pl.show()

