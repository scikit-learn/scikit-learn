"""
=========================================
Plot PCA 2d projection of of Iris dataset
=========================================
"""

import pylab as pl

from scikits.learn import datasets
from scikits.learn.pca import PCA

iris = datasets.load_iris()

X = iris.data
y = iris.target
target_names = iris.target_names

pca = PCA(k=2)
X_r = pca.fit(X).transform(X)

# Percentage of variance explained for each components
print pca.explained_variance_

pl.figure()
pl.hold('on')
for c, i, target_name in zip("rgb", [0, 1, 2], target_names):
   pl.scatter(X_r[y==i,0], X_r[y==i,1], c=c, label=target_name)
pl.hold('off')
pl.legend()
pl.title('PCA of IRIS dataset')

pl.show()

