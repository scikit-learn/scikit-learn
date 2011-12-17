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

print y_pred
print y

