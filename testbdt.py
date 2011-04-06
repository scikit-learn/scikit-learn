from scikits.learn.decisiontree import DecisionTree
from scikits.learn.ensemble.boosting import AdaBoost
import numpy as np
import matplotlib.pyplot as plt

#a = AdaBoost(DecisionTree, minleafsize=1000, nbins=1000)
a = DecisionTree(minleafsize=100)

cov = [[1,0],[0,1]]
X_bkg = np.reshape(np.random.multivariate_normal([-1,-1], cov, 100000), (100000,-1))
X_sig = np.reshape(np.random.multivariate_normal([1,1], cov, 100000), (100000, -1))

x = np.concatenate([X_sig, X_bkg])
y = np.append(np.ones(100000),-2*np.ones(100000))

a.fit(x, y)

plt.figure()
plt.hist(a.predict(x))

plt.figure()
Xs,Ys = X_sig.T
Xb,Yb = X_bkg.T
plt.plot(Xs,Ys,'x')
plt.plot(Xb,Yb,'o',alpha=.3)
plt.axis('equal')
plt.show()
