"""
================================================================
Plot the decision tree
================================================================

Plot the decision tree trained on some features of the iris dataset.

See :ref:`decision tree <tree>` for more information on the estimator.
"""
print(__doc__)

from sklearn import tree
from sklearn.datasets import load_iris

# Load data
iris = load_iris()

# We only take the two parameters.
X, y = iris.data[:, 2:4], iris.target

# Train
clf = tree.DecisionTreeClassifier().fit(X, y)

# And plot
tree.plot(clf, height=700, width=500, filled=True, rounded=True,
          special_characters=True, proportion=True)
