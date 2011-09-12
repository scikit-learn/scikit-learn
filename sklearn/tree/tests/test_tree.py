"""
Testing for Tree module (sklearn.tree)

"""

import numpy as np
from numpy.testing import assert_array_equal, assert_array_almost_equal, \
                          assert_almost_equal
from nose.tools import assert_raises, with_setup

from sklearn import tree, datasets, metrics

# toy sample
X = [[-2, -1], [-1, -1], [-1, -2], [1, 1], [1, 2], [2, 1]]
Y = [-1, -1, -1, 1, 1, 1]
T = [[-1, -1], [2, 2], [3, 2]]
true_result = [-1, 1, 1]

TEMP_DIR = ""

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

    clf = tree.DecisionTreeClassifier()
    clf.fit(X, Y)

    assert_array_equal(clf.predict(T), true_result)

    # With subsampling
    clf = tree.DecisionTreeClassifier(max_features=1, random_state=1)
    clf.fit(X, Y)

    assert_array_equal(clf.predict(T), true_result)

    # With n_classes given
    clf = tree.DecisionTreeClassifier(n_classes=2)
    clf.fit(X, Y)

    assert_array_equal(clf.predict(T), true_result)


def test_regression_toy():
    """Check regression on a toy dataset."""
    clf = tree.DecisionTreeRegressor()
    clf.fit(X, Y)

    assert_almost_equal(clf.predict(T), true_result)

    # With subsampling
    clf = tree.DecisionTreeRegressor(max_features=1, random_state=1)
    clf.fit(X, Y)

    assert_almost_equal(clf.predict(T), true_result)


def setup_graphviz_toy():
    import tempfile
    global TEMP_DIR
    TEMP_DIR = tempfile.mkdtemp(dir=".")


def teardown_graphviz_toy():
    import shutil
    shutil.rmtree(TEMP_DIR, ignore_errors=True)


@with_setup(setup_graphviz_toy, teardown_graphviz_toy)
def test_graphviz_toy():
    """Check correctness of graphviz output on a toy dataset."""
    clf = tree.DecisionTreeClassifier(max_depth=100, min_split=1)
    clf.fit(X, Y)
    out = open(TEMP_DIR + "/tree.dot", 'w')
    exporter = tree.GraphvizExporter(out)
    clf.export(exporter)
    exporter.close()

    import os
    dirname = os.path.dirname(__file__)
    if dirname != "":
        dirname += "/"
    with open(TEMP_DIR + "/tree.dot") as f1:
        with open(dirname + 'test_tree.dot') as f2:
            # replace unique memory addresses with a tmp string
            l1 = f1.read()
            l2 = f2.read()
            assert l1 == l2, \
                "graphviz output test failed\n: %s != %s" % (l1, l2)


def test_iris():
    """Check consistency on dataset iris."""

    for c in ('gini', \
              'entropy'):
        clf = tree.DecisionTreeClassifier(criterion=c)\
              .fit(iris.data, iris.target)

        score = np.mean(clf.predict(iris.data) == iris.target)
        assert score > 0.9, "Failed with criterion " + c + \
            " and score = " + str(score)

        clf = tree.DecisionTreeClassifier(criterion=c,
                                          max_features=2,
                                          random_state=1)\
              .fit(iris.data, iris.target)

        score = np.mean(clf.predict(iris.data) == iris.target)
        assert score > 0.5, "Failed with criterion " + c + \
            " and score = " + str(score)


def test_boston():
    """Check consistency on dataset boston house prices."""
    for c in ('mse',):
        clf = tree.DecisionTreeRegressor(criterion=c)\
              .fit(boston.data, boston.target)

        score = np.mean(np.power(clf.predict(boston.data) - boston.target, 2))
        assert score < 1, "Failed with criterion " + c + \
            " and score = " + str(score)

        clf = tree.DecisionTreeRegressor(criterion=c,
                                         max_features=6,
                                         random_state=1)\
              .fit(boston.data, boston.target)

        #using fewer features reduces the learning ability of this tree,
        # but reduces training time.
        score = np.mean(np.power(clf.predict(boston.data) - boston.target, 2))
        assert score < 2, "Failed with criterion " + c + \
            " and score = " + str(score)


def test_probability():
    """Predict probabilities using DecisionTreeClassifier."""

    clf = tree.DecisionTreeClassifier()
    clf.fit(iris.data, iris.target)

    prob_predict = clf.predict_proba(iris.data)
    assert_array_almost_equal(
        np.sum(prob_predict, 1), np.ones(iris.data.shape[0]))
    assert np.mean(np.argmax(prob_predict, 1)
                   == clf.predict(iris.data)) > 0.9

    assert_almost_equal(clf.predict_proba(iris.data),
                        np.exp(clf.predict_log_proba(iris.data)), 8)


def test_error():
    """Test that it gives proper exception on deficient input."""
    # impossible value of min_split
    assert_raises(ValueError, \
                  tree.DecisionTreeClassifier(min_split=-1).fit, X, Y)

    # impossible value of max_depth
    assert_raises(ValueError, \
                  tree.DecisionTreeClassifier(max_depth=-1).fit, X, Y)

    clf = tree.DecisionTreeClassifier()

    Y2 = Y[:-1]  # wrong dimensions for labels
    assert_raises(ValueError, clf.fit, X, Y2)

    # Test with arrays that are non-contiguous.
    Xf = np.asfortranarray(X)
    clf = tree.DecisionTreeClassifier()
    clf.fit(Xf, Y)
    assert_array_equal(clf.predict(T), true_result)

    # use values of max_features that are invalid
    clf = tree.DecisionTreeClassifier(max_features=-1)
    assert_raises(ValueError, clf.fit, X, Y2)

    clf = tree.DecisionTreeClassifier(max_features=10)
    assert_raises(ValueError, clf.fit, X, Y2)

    clf = tree.DecisionTreeClassifier()
    # predict before fitting
    assert_raises(Exception, clf.predict, T)

    # predict on vector with different dims
    clf.fit(X, Y)
    t = np.asanyarray(T)
    assert_raises(ValueError, clf.predict, t[:, 1:])

    # labels out of range
    clf = tree.DecisionTreeClassifier(n_classes=1)
    assert_raises(ValueError, clf.fit, X, Y2)

    clf = tree.DecisionTreeClassifier(n_classes=3)
    assert_raises(ValueError, clf.fit, X, Y)

    # max_features invalid
    clf = tree.DecisionTreeClassifier(max_features=-1)
    assert_raises(ValueError, clf.fit, X, Y)

    clf = tree.DecisionTreeClassifier(max_features=3)
    assert_raises(ValueError, clf.fit, X, Y)

    # predict before fit
    clf = tree.DecisionTreeClassifier()
    assert_raises(Exception, clf.predict_proba, X)

    clf.fit(X, Y)
    X2 = [-2, -1, 1]  # wrong feature shape for sample
    assert_raises(ValueError, clf.predict_proba, X2)

    Xt = np.array(X).T

    clf = tree.DecisionTreeClassifier()
    clf.fit(np.dot(X, Xt), Y)
    assert_raises(ValueError, clf.predict, X)

    clf = tree.DecisionTreeClassifier()
    clf.fit(X, Y)
    assert_raises(ValueError, clf.predict, Xt)


if __name__ == '__main__':
    import nose
    nose.runmodule()
