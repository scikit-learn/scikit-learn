"""
Simple tests to demonstrate usage.
Proper test module will come.
"""

from scikits.learn.decisiontree import DecisionTree
from scikits.learn.ensemble.boosting import AdaBoost
from scikits.learn.ensemble.bagging import Bagged
import numpy as np
import matplotlib.pyplot as plt

boosted = AdaBoost(DecisionTree, minleafsize=100, nbins = 100)
bagged = Bagged(DecisionTree, minleafsize=100, nbins = 100)
# use maxdepth=1 for decision stumps

nsample = 10000
#cov = [[.5,0,0],[0,.5,0],[0,0,0.5]]
cov = [[.5,0],[0,.5]]
X_bkg = np.concatenate([np.reshape(np.random.multivariate_normal([0,2], cov, nsample/2), (nsample/2,-1)),
                        np.reshape(np.random.multivariate_normal([0,-2], cov, nsample/2), (nsample/2,-1))])
X_sig = np.reshape(np.random.multivariate_normal([-1.5,0], cov, nsample), (nsample,-1))

x = np.concatenate([X_sig, X_bkg])
y = np.append(np.ones(nsample),-np.ones(nsample))

boosted.fit(x, y, boosts = 20, beta = .2)
bagged.fit(x, y, baggs = 20, sample_fraction = .6)

plt.figure()
plt.hist(boosted.predict(X_bkg), bins=20, range=(-1,1), label="Boosted Background")
plt.hist(boosted.predict(X_sig), bins=20, range=(-1,1), label="Boosted Signal", alpha=.5)
l = plt.legend()
l.legendPatch.set_alpha(0.5)

plt.figure()
plt.hist(bagged.predict(X_bkg), bins=20, range=(-1,1), label="Bagged Background")
plt.hist(bagged.predict(X_sig), bins=20, range=(-1,1), label="Bagged Signal", alpha=.5)
l = plt.legend()
l.legendPatch.set_alpha(0.5)

plt.figure()
Xs,Ys = X_sig.T
Xb,Yb = X_bkg.T
plt.plot(Xb,Yb,'o', label="Background")
plt.plot(Xs,Ys,'o', label="Signal", alpha=.5)
l =plt.legend()
l.legendPatch.set_alpha(0.5)
plt.axis('equal')
plt.show()
