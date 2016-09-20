"""
================================================================
Plot the decision tree
================================================================

Plot the decision tree trained on some features of the iris dataset.

See :ref:`decision tree <tree>` for more information on the estimator.
"""
print(__doc__)

try:
    import graphviz
except ImportError:
    raise ImportError('This example requires graphviz to be installed')

import six

from matplotlib import pyplot
from scipy.misc import imread
from sklearn import tree
from sklearn.datasets import load_iris

# Load data
iris = load_iris()

# We only take the two parameters.
X, y = iris.data[:, 2:4], iris.target

# Train
clf = tree.DecisionTreeClassifier().fit(X, y)

# And plot
tree_source = graphviz.Source(tree.get_graphviz_source(
    clf, filled=True, rounded=True, special_characters=True, proportion=True))
image = imread(six.BytesIO(tree_source.pipe(format='png')))
pyplot.imshow(image)
pyplot.show()
