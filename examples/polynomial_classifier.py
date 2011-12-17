"""

Polynomial Classifier - An example

"""

# data
from sklearn import datasets
iris = datasets.load_iris()

X = iris.data
y = iris.target

# apply classifier
from sklearn.polynomial_classifier import PC
pc = PC(degree=2)

y_pred = pc.fit(X, y).predict(X)

print "Number of mislabeled points : %d" % (y != y_pred).sum()
