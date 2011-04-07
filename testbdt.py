from scikits.learn.decisiontree import DecisionTree
from scikits.learn.ensemble.boosting import AdaBoost
import numpy as np
import matplotlib.pyplot as plt

a = AdaBoost(DecisionTree, minleafsize=100, maxdepth = 1, nbins=100)
#a = DecisionTree(minleafsize=100, maxdepth=1, nbins=100)

nsample = 50000
#cov = [[.5,0,0],[0,.5,0],[0,0,0.5]]
cov = [[.5,0],[0,.5]]
X_bkg = np.reshape(np.random.multivariate_normal([.5,.5], cov, nsample), (nsample,-1))
X_sig = np.reshape(np.random.multivariate_normal([-.5,-.5], cov, nsample), (nsample,-1))

x = np.concatenate([X_sig, X_bkg])
y = np.append(np.ones(nsample),-np.ones(nsample))

a.fit(x, y, boosts=20)

plt.figure()
plt.hist(a.predict(X_bkg), bins=20, range=(-1,1), label="Background")
plt.hist(a.predict(X_sig), bins=20, range=(-1,1), label="Signal", alpha=.5)
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
