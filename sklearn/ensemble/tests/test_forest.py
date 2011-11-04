import numpy as np
from numpy.testing import assert_array_equal, assert_array_almost_equal, assert_almost_equal

from sklearn import datasets, ensemble

# toy sample
X = [[-2, -1], [-1, -1], [-1, -2], [1, 1], [1, 2], [2, 1]]
y = [-1, -1, -1, 1, 1, 1]
T = [[-1, -1], [2, 2], [3, 2]]
true_result = [-1, 1, 1]

# also load the iris dataset
# and randomly permute it
iris = datasets.load_iris()
np.random.seed([1])
perm = np.random.permutation(iris.target.size)
iris.data = iris.data[perm]
iris.target = iris.target[perm]

# also load the boston dataset
# and randomly permute it
boston = datasets.load_boston()
perm = np.random.permutation(boston.target.size)
boston.data = boston.data[perm]
boston.target = boston.target[perm]

def test_classification_toy():
    """Check classification on a toy dataset."""
    clf = ensemble.RandomForestClassifier(n_trees=5)
    clf.fit(X, y)

    assert_array_equal(clf.predict(T), true_result)

    clf = ensemble.RandomForestClassifier(n_trees=5, max_features=1)
    clf.fit(X, y)

    assert_array_equal(clf.predict(T), true_result)
