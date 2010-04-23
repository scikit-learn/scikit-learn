"""
==================
One-class SVM
==================
"""

import numpy as np
import pylab as pl
from scikits.learn import svm

xx, yy = np.meshgrid(np.linspace(-5, 5, 500), np.linspace(-5, 5, 500))
X = np.random.randn(100, 2)
Y = [0]*100

# fit the model
clf = svm.OneClassSVM(nu=0.5)
clf.fit(X, Y)

# plot the line, the points, and the nearest vectors to the plane
Z = clf.predict(np.c_[xx.ravel(), yy.ravel()])
Z = Z.reshape(xx.shape)

pl.set_cmap(pl.cm.Paired)
pl.pcolormesh(xx, yy, Z)
pl.scatter(X[:,0], X[:,1], c=Y)
pl.scatter(clf.support_[:,0], clf.support_[:,1], c='black')
pl.axis('tight')
pl.show()
