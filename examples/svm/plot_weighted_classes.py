"""
================================================
SVM: Separating hyperplane with weighted classes
================================================

"""

import numpy as np
import pylab as pl
from scikits.learn import svm

# we create 40 separable points
np.random.seed(0)
n_samples_1 = 1000
n_samples_2 = 100
X = np.r_[1.5*np.random.randn(n_samples_1, 2), 0.5*np.random.randn(n_samples_2, 2) + [2, 2]]
Y = [0]*(n_samples_1) + [1]*(n_samples_2)

# fit the model and get the separating hyperplane
clf = svm.SVC(kernel='linear')
clf.fit(X, Y)

w = clf.coef_[0]
a = -w[0]/w[1]
xx = np.linspace(-5, 5)
yy = a*xx - (clf.intercept_[0])/w[1]


# get the separating hyperplane using weighted classes
wclf = svm.SVC(kernel='linear')
wclf.fit(X, Y, {1: 10})

ww = wclf.coef_[0]
wa = -ww[0]/ww[1]
wyy = wa*xx - (wclf.intercept_[0])/ww[1]

# plot separating hyperplanes and samples
pl.set_cmap(pl.cm.Paired)
pl.plot(xx, yy, 'k-')
pl.plot(xx, wyy, 'k--')
pl.scatter(X[:,0], X[:,1], c=Y)

pl.axis('tight')
pl.show()

