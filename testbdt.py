from scikits.learn.decisiontree import DecisionTree
from scikits.learn.ensemble.boosting import AdaBoost
import numpy as np
import matplotlib.pyplot as plt

a = AdaBoost(DecisionTree, minleafsize=2, maxdepth = 1)
#a = DecisionTree(minleafsize=100, maxdepth=1)

cov = [[.5,0],[0,.5]]
X_bkg = np.reshape(np.random.multivariate_normal([-1,-1], cov, 10), (10,-1))
X_sig = np.reshape(np.random.multivariate_normal([1,1], cov, 10), (10, -1))

x = np.concatenate([X_sig, X_bkg])
y = np.append(np.ones(10),-np.ones(10))

a.fit(x, y, boosts = 10)

plt.figure()
plt.hist(a.predict(x))

plt.figure()
Xs,Ys = X_sig.T
Xb,Yb = X_bkg.T
plt.plot(Xs,Ys,'x')
plt.plot(Xb,Yb,'o',alpha=.3)
plt.axis('equal')
#plt.show()
