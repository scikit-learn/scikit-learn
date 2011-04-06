from scikits.learn.decisiontree import DecisionTree
import numpy as np

a = DecisionTree(minleafsize=1)

x = np.array([[.3, .4, .2] for i in xrange(100)])
print x.shape
y = np.append(np.ones(50),-np.ones(50))

a.fit(x,y)
