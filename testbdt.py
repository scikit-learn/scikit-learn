from scikits.learn.decisiontree import DecisionTree
from scikits.learn.ensemble.boosting import AdaBoost
import numpy as np
import matplotlib.pyplot as plt

a = AdaBoost(DecisionTree, minleafsize=100, maxdepth = 1)
#a = DecisionTree(minleafsize=100, maxdepth=1)

nsample = 10000
cov = [[.5,0,0],[0,.5,0],[0,0,0.5]]
X_bkg = np.reshape(np.random.multivariate_normal([-2,-2,0], cov, nsample), (nsample,-1))
X_sig = np.reshape(np.random.multivariate_normal([2,2,0], cov, nsample), (nsample,-1))


x = np.concatenate([X_sig, X_bkg])
y = np.append(np.ones(nsample),-np.ones(nsample))

a.fit(x, y, boosts = 100)

plt.figure()
p = a.predict(x)
plt.hist(p, bins=100)

plt.figure()
Xs,Ys,Zs = X_sig.T
Xb,Yb,Zb = X_bkg.T
plt.plot(Xs,Ys,'x')
plt.plot(Xb,Yb,'o',alpha=.3)
plt.axis('equal')
plt.show()
