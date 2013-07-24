"""
==============================================
Plot randomly generated classification dataset
==============================================

Plot several randomly generated 2D classification datasets.
This example illustrates the `datasets.make_classification`
function.

Three binary and two multi-class classification datasets
are generated, with different numbers of informative
features and clusters per class.
"""

print(__doc__)

import pylab as pl

from sklearn.datasets import make_classification

pl.figure(figsize=(8, 6))
pl.subplots_adjust(bottom=.05, top=.9, left=.05, right=.95)

pl.subplot(221)
pl.title("One informative feature, one cluster", fontsize='small')
X1, Y1 = make_classification(n_features=2, n_redundant=0, n_informative=1,
                             n_clusters_per_class=1)
pl.scatter(X1[:, 0], X1[:, 1], marker='o', c=Y1)

pl.subplot(222)
pl.title("Two informative features, one cluster", fontsize='small')
X1, Y1 = make_classification(n_features=2, n_redundant=0, n_informative=2,
                             n_clusters_per_class=1)
pl.scatter(X1[:, 0], X1[:, 1], marker='o', c=Y1)

pl.subplot(223)
pl.title("Two informative features, two clusters", fontsize='small')
X2, Y2 = make_classification(n_features=2, n_redundant=0, n_informative=2)
pl.scatter(X2[:, 0], X2[:, 1], marker='o', c=Y2)


pl.subplot(224)
pl.title("Multi-class, two informative features, one cluster",
         fontsize='small')
X1, Y1 = make_classification(n_features=2, n_redundant=0, n_informative=2,
                             n_clusters_per_class=1, n_classes=3)
pl.scatter(X1[:, 0], X1[:, 1], marker='o', c=Y1)

pl.show()
