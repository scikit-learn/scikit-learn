"""
Testing for the boost module (sklearn.ensemble.boost).
"""

import numpy as np
from numpy.testing import assert_array_equal
from numpy.testing import assert_array_almost_equal
from numpy.testing import assert_equal
from nose.tools import assert_true
from nose.tools import assert_raises

from sklearn.grid_search import GridSearchCV
from sklearn.ensemble import AdaBoostClassifier, AdaBoostRegressor
from sklearn.tree import DecisionTreeClassifier
from sklearn import datasets

# toy sample
X = [[-2, -1], [-1, -1], [-1, -2], [1, 1], [1, 2], [2, 1]]
y = [-1, -1, -1, 1, 1, 1]
T = [[-1, -1], [2, 2], [3, 2]]
true_result = [-1, 1, 1]

# also load the iris dataset
# and randomly permute it
iris = datasets.load_iris()
rng = np.random.RandomState(0)
perm = rng.permutation(iris.target.size)
iris.data = iris.data[perm]
iris.target = iris.target[perm]

# also load the boston dataset
# and randomly permute it
boston = datasets.load_boston()
perm = rng.permutation(boston.target.size)
boston.data = boston.data[perm]
boston.target = boston.target[perm]


def test_classification_toy():
    """Check classification on a toy dataset."""
    clf = AdaBoostClassifier(n_estimators=10)
    clf.fit(X, y)
    assert_array_equal(clf.predict(T), true_result)


def test_regression_toy():
    """Check classification on a toy dataset."""
    clf = AdaBoostRegressor(n_estimators=10)
    clf.fit(X, y)
    assert_array_equal(clf.predict(T), true_result)


def test_iris():
    """Check consistency on dataset iris."""
    for c in ("gini", "entropy"):
        # AdaBoost classification
        clf = AdaBoostClassifier(DecisionTreeClassifier(criterion=c),
                                 n_estimators=1)
        clf.fit(iris.data, iris.target)
        score = clf.score(iris.data, iris.target)
        assert score > 0.9, "Failed with criterion %s and score = %f" % (c,
                                                                         score)


def test_boston():
    """Check consistency on dataset boston house prices."""
    clf = AdaBoostRegressor(n_estimators=50)
    clf.fit(boston.data, boston.target)
    score = clf.score(boston.data, boston.target)
    assert score > 0.85


def test_probability():
    """Predict probabilities."""
    # AdaBoost classification
    clf = AdaBoostClassifier(n_estimators=10)
    clf.fit(iris.data, iris.target)

    assert_array_almost_equal(np.sum(clf.predict_proba(iris.data), axis=1),
                              np.ones(iris.data.shape[0]))
    assert_array_almost_equal(clf.predict_proba(iris.data),
                              np.exp(clf.predict_log_proba(iris.data)))


def test_staged_predict():
    """Check staged predictions."""
    # AdaBoost classification
    clf = AdaBoostClassifier(n_estimators=10)
    clf.fit(iris.data, iris.target)

    predictions = clf.predict(iris.data)
    staged_predictions = [p for p in clf.staged_predict(iris.data)]
    proba = clf.predict_proba(iris.data)
    staged_probas = [p for p in clf.staged_predict_proba(iris.data)]
    score = clf.score(iris.data, iris.target)
    staged_scores = [s for s in clf.staged_score(iris.data, iris.target)]

    assert_equal(len(staged_predictions), 10)
    assert_array_equal(predictions, staged_predictions[-1])
    assert_equal(len(staged_probas), 10)
    assert_array_equal(proba, staged_probas[-1])
    assert_equal(len(staged_scores), 10)
    assert_array_equal(score, staged_scores[-1])

    # AdaBoost regression
    clf = AdaBoostRegressor(n_estimators=10)
    clf.fit(boston.data, boston.target)

    predictions = clf.predict(boston.data)
    staged_predictions = [p for p in clf.staged_predict(boston.data)]
    score = clf.score(boston.data, boston.target)
    staged_scores = [s for s in clf.staged_score(boston.data, boston.target)]

    assert_equal(len(staged_predictions), 10)
    assert_array_equal(predictions, staged_predictions[-1])
    assert_equal(len(staged_scores), 10)
    assert_array_equal(score, staged_scores[-1])


def test_gridsearch():
    """Check that base trees can be grid-searched."""
    # AdaBoost classification
    boost = AdaBoostClassifier()
    parameters = {'n_estimators': (1, 2),
                  'base_estimator__max_depth': (1, 2)}
    clf = GridSearchCV(boost, parameters)
    clf.fit(iris.data, iris.target)

    # AdaBoost regression
    boost = AdaBoostRegressor()
    parameters = {'n_estimators': (1, 2),
                  'base_estimator__max_depth': (1, 2)}
    clf = GridSearchCV(boost, parameters)
    clf.fit(boston.data, boston.target)


def test_pickle():
    """Check pickability."""
    import pickle

    # Adaboost classifier
    obj = AdaBoostClassifier()
    obj.fit(iris.data, iris.target)
    score = obj.score(iris.data, iris.target)
    s = pickle.dumps(obj)

    obj2 = pickle.loads(s)
    assert_equal(type(obj2), obj.__class__)
    score2 = obj2.score(iris.data, iris.target)
    assert score == score2

    # Adaboost regressor
    obj = AdaBoostRegressor()
    obj.fit(boston.data, boston.target)
    score = obj.score(boston.data, boston.target)
    s = pickle.dumps(obj)

    obj2 = pickle.loads(s)
    assert_equal(type(obj2), obj.__class__)
    score2 = obj2.score(boston.data, boston.target)
    assert score == score2


def test_importances():
    """Check variable importances."""
    X, y = datasets.make_classification(n_samples=2000,
                                        n_features=10,
                                        n_informative=3,
                                        n_redundant=0,
                                        n_repeated=0,
                                        shuffle=False,
                                        random_state=1)

    clf = AdaBoostClassifier(compute_importances=True, n_estimators=50)
    clf.fit(X, y)
    importances = clf.feature_importances_
    n_important = sum(importances > 0.1)

    assert_equal(importances.shape[0], 10)
    assert_equal(n_important, 3)

    clf = AdaBoostClassifier()
    clf.fit(X, y)
    assert_true(clf.feature_importances_ is None)


def test_error():
    """Test that it gives proper exception on deficient input."""
    # Invalid values for parameters
    assert_raises(ValueError,
                  AdaBoostClassifier(learning_rate=-1).fit,
                  X, y)


def test_base_estimator():
    """Test different base estimators."""
    from sklearn.tree import DecisionTreeClassifier
    from sklearn.ensemble import RandomForestClassifier
    from sklearn.svm import SVC

    clf = AdaBoostClassifier(base_estimator=DecisionTreeClassifier())
    clf.fit(X, y)

    clf = AdaBoostClassifier(base_estimator=RandomForestClassifier())
    clf.fit(X, y)

    clf = AdaBoostClassifier(base_estimator=SVC())
    clf.fit(X, y)

    from sklearn.tree import DecisionTreeRegressor
    from sklearn.ensemble import RandomForestRegressor
    from sklearn.svm import SVR

    clf = AdaBoostRegressor(base_estimator=DecisionTreeRegressor())
    clf.fit(X, y)

    clf = AdaBoostRegressor(base_estimator=RandomForestRegressor())
    clf.fit(X, y)

    clf = AdaBoostRegressor(base_estimator=SVR())
    clf.fit(X, y)


if __name__ == "__main__":
    import nose
    nose.runmodule()
