# Author: Alexandre Gramfort <alexandre.gramfort@inria.fr>
# License: BSD Style.

# $Id$

import pylab as pl
import numpy as np

from scikits.learn.logistic import LogisticRegression
from scikits.learn.svm import SVC
from scikits.learn import datasets

iris = datasets.load_iris()
X = iris.data[:, :2] # we only take the first two features for visualization
y = iris.target

n_features = X.shape[1]

C = 1.0

# Create classifier (any of the following 3)
classifier = LogisticRegression(C=C, penalty='l1')
classifier = LogisticRegression(C=C, penalty='l2')
classifier = SVC(kernel='linear', C=C, probability=True)

classifier.fit(X, y)

y_pred = classifier.predict(X)
classif_rate = np.mean(y_pred.ravel() == y.ravel()) * 100
print  "classif_rate : %f " % classif_rate

# ======================
# = View probabilities =
# ======================
pl.figure()
xx = np.linspace(3,9,100)
yy = np.linspace(1,5,100).T
xx, yy = np.meshgrid(xx, yy)
Xfull = np.c_[xx.ravel(),yy.ravel()]
probas = classifier.predict_proba(Xfull)
n_classes = np.unique(y_pred).size
for k in range(n_classes):
    pl.subplot(1, n_classes, k + 1)
    pl.title("Class %d" % k)
    imshow_handle = pl.imshow(probas[:,k].reshape((100, 100)), extent=(3, 9, 1, 5), origin='lower')
    pl.hold(True)
    idx = (y_pred == k)
    if idx.any(): pl.scatter(X[idx,0], X[idx,1], marker='o', c='k')

ax = pl.axes([0.15,0.04,0.7,0.05])
pl.title("Probability")
pl.colorbar(imshow_handle, cax=ax, orientation='horizontal')

pl.show()
