import numpy as np

from sklearn.utils.testing import assert_array_equal

from sklearn.random_classifier import RandomClassifier


def test_most_frequent_sampling():
    X = [[0], [0], [0], [0]] # ignored
    y = [1, 2, 1, 1]

    clf = RandomClassifier(sampling="most_frequent", random_state=0)
    clf.fit(X, y)
    assert_array_equal(clf.predict(X, y), np.ones(len(y)))

def test_stratified_sampling():
    X = [[0], [0], [0], [0], [0]] # ignored
    y = [1, 2, 1, 1, 2]

    clf = RandomClassifier(sampling="stratified", random_state=0)
    clf.fit(X, y)
    assert_array_equal(clf.predict(X, y), [1, 2, 2, 1, 1])


def test_uniform_sampling():
    X = [[0], [0], [0], [0]] # ignored
    y = [1, 2, 1, 1]

    clf = RandomClassifier(sampling="uniform", random_state=0)
    clf.fit(X, y)
    assert_array_equal(clf.predict(X, y), [1, 2, 2, 1])

