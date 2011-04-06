from scikits.learn.decisiontree import DecisionTree
import numpy as np

a = DecisionTree(minleafsize=10)

cov = [[1,.5],[.5,1]]
X_bkg = np.random.multivariate_normal([.5,.5], cov, 1000)
X_sig = np.random.multivariate_normal([-.5,-.5], cov, 1000)

x = np.append(X_sig, X_bkg)
x = np.reshape(x, (2000,-1))
print x.shape
y = np.append(np.ones(1000),-np.ones(1000))
print y.shape

a.fit(x,y)

print a.predict(x)
