import numpy as np
import pylab as pl
from scikits.learn import svm

# we create 40 separable points
np.random.seed(0)
X = np.r_[np.random.randn(20, 2) - [2,2], np.random.randn(20, 2) + [2, 2]]
Y = [0]*20 + [1]*20

# fit the model
clf = svm.SVC(kernel='linear')
clf.fit(X, Y)

# get the separating hyperplane
w = np.dot(clf.coef_[0], clf.support_) 
xx = np.linspace(-2, 2)
yy = (clf.rho_[0] - w[0]*xx)/w[1]

# plot the line, the points, and the nearest vectors to the plane
pl.set_cmap(pl.cm.Paired)
pl.plot(xx, yy, 'k--')
pl.scatter(X[:,0], X[:,1], c=Y)
pl.scatter(clf.support_[:,0], clf.support_[:,1], marker='+')

pl.axis('tight')
pl.show()

