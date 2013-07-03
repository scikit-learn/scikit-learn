"""
Testing for the boost module (sklearn.ensemble.boost).
"""

import numpy as np
from numpy.testing import assert_array_equal, assert_array_less
from numpy.testing import assert_array_almost_equal
from numpy.testing import assert_equal
from nose.tools import assert_raises

from sklearn.dummy import DummyClassifier
from sklearn.dummy import DummyRegressor
from sklearn.grid_search import GridSearchCV
from sklearn.ensemble import AdaBoostClassifier
from sklearn.ensemble import AdaBoostRegressor
from sklearn.utils import shuffle
from sklearn import datasets


# Common random state
rng = np.random.RandomState(0)

# Toy sample
X = [[-2, -1], [-1, -1], [-1, -2], [1, 1], [1, 2], [2, 1]]
y_class = ["foo", "foo", "foo", 1, 1, 1]    # test string class labels
y_regr = [-1, -1, -1, 1, 1, 1]
T = [[-1, -1], [2, 2], [3, 2]]
y_t_class = ["foo", 1, 1]
y_t_regr = [-1, 1, 1]

# Load the iris dataset and randomly permute it
iris = datasets.load_iris()
perm = rng.permutation(iris.target.size)
iris.data, iris.target = shuffle(iris.data, iris.target, random_state=rng)

# Load the boston dataset and randomly permute it
boston = datasets.load_boston()
boston.data, boston.target = shuffle(boston.data, boston.target,
                                     random_state=rng)


def test_classification_toy():
    """Check classification on a toy dataset."""
    for alg in ['SAMME', 'SAMME.R']:
        clf = AdaBoostClassifier(algorithm=alg)
        clf.fit(X, y_class)
        assert_array_equal(clf.predict(T), y_t_class)
        assert_array_equal(np.unique(np.asarray(y_t_class)), clf.classes_)
        assert_equal(clf.predict_proba(T).shape, (len(T), 2))
        assert_equal(clf.decision_function(T).shape, (len(T),))


def test_regression_toy():
    """Check classification on a toy dataset."""
    clf = AdaBoostRegressor(random_state=0)
    clf = AdaBoostRegressor()
    clf.fit(X, y_regr)
    assert_array_equal(clf.predict(T), y_t_regr)


def test_iris():
    """Check consistency on dataset iris."""
    classes = np.unique(iris.target)
    clf_samme = prob_samme = None

    for alg in ['SAMME', 'SAMME.R']:
        clf = AdaBoostClassifier(algorithm=alg)
        clf.fit(iris.data, iris.target)

        assert_array_equal(classes, clf.classes_)
        proba = clf.predict_proba(iris.data)
        if alg == "SAMME":
            clf_samme = clf
            prob_samme = proba
        assert_equal(proba.shape[1], len(classes))
        assert_equal(clf.decision_function(iris.data).shape[1], len(classes))

        score = clf.score(iris.data, iris.target)
        assert score > 0.9, "Failed with algorithm %s and score = %f" % \
            (alg, score)

    # Somewhat hacky regression test: prior to
    # ae7adc880d624615a34bafdb1d75ef67051b8200,
    # predict_proba returned SAMME.R values for SAMME.
    clf_samme.algorithm = "SAMME.R"
    assert_array_less(0,
                      np.abs(clf_samme.predict_proba(iris.data) - prob_samme))


def test_boston():
    """Check consistency on dataset boston house prices."""
    clf = AdaBoostRegressor(random_state=0)
    clf.fit(boston.data, boston.target)
    score = clf.score(boston.data, boston.target)
    assert score > 0.85


def test_staged_predict():
    """Check staged predictions."""
    # AdaBoost classification
    for alg in ['SAMME', 'SAMME.R']:
        clf = AdaBoostClassifier(algorithm=alg, n_estimators=10)
        clf.fit(iris.data, iris.target)

        predictions = clf.predict(iris.data)
        staged_predictions = [p for p in clf.staged_predict(iris.data)]
        proba = clf.predict_proba(iris.data)
        staged_probas = [p for p in clf.staged_predict_proba(iris.data)]
        score = clf.score(iris.data, iris.target)
        staged_scores = [s for s in clf.staged_score(iris.data, iris.target)]

        assert_equal(len(staged_predictions), 10)
        assert_array_almost_equal(predictions, staged_predictions[-1])
        assert_equal(len(staged_probas), 10)
        assert_array_almost_equal(proba, staged_probas[-1])
        assert_equal(len(staged_scores), 10)
        assert_array_almost_equal(score, staged_scores[-1])

    # AdaBoost regression
    clf = AdaBoostRegressor(n_estimators=10, random_state=0)
    clf.fit(boston.data, boston.target)

    predictions = clf.predict(boston.data)
    staged_predictions = [p for p in clf.staged_predict(boston.data)]
    score = clf.score(boston.data, boston.target)
    staged_scores = [s for s in clf.staged_score(boston.data, boston.target)]

    assert_equal(len(staged_predictions), 10)
    assert_array_almost_equal(predictions, staged_predictions[-1])
    assert_equal(len(staged_scores), 10)
    assert_array_almost_equal(score, staged_scores[-1])


def test_gridsearch():
    """Check that base trees can be grid-searched."""
    # AdaBoost classification
    boost = AdaBoostClassifier()
    parameters = {'n_estimators': (1, 2),
                  'base_estimator__max_depth': (1, 2),
                  'algorithm': ('SAMME', 'SAMME.R')}
    clf = GridSearchCV(boost, parameters)
    clf.fit(iris.data, iris.target)

    # AdaBoost regression
    boost = AdaBoostRegressor(random_state=0)
    parameters = {'n_estimators': (1, 2),
                  'base_estimator__max_depth': (1, 2)}
    clf = GridSearchCV(boost, parameters)
    clf.fit(boston.data, boston.target)


def test_pickle():
    """Check pickability."""
    import pickle

    # Adaboost classifier
    for alg in ['SAMME', 'SAMME.R']:
        obj = AdaBoostClassifier(algorithm=alg)
        obj.fit(iris.data, iris.target)
        score = obj.score(iris.data, iris.target)
        s = pickle.dumps(obj)

        obj2 = pickle.loads(s)
        assert_equal(type(obj2), obj.__class__)
        score2 = obj2.score(iris.data, iris.target)
        assert_equal(score, score2)

    # Adaboost regressor
    obj = AdaBoostRegressor(random_state=0)
    obj.fit(boston.data, boston.target)
    score = obj.score(boston.data, boston.target)
    s = pickle.dumps(obj)

    obj2 = pickle.loads(s)
    assert_equal(type(obj2), obj.__class__)
    score2 = obj2.score(boston.data, boston.target)
    assert_equal(score, score2)


def test_importances():
    """Check variable importances."""
    X, y = datasets.make_classification(n_samples=2000,
                                        n_features=10,
                                        n_informative=3,
                                        n_redundant=0,
                                        n_repeated=0,
                                        shuffle=False,
                                        random_state=1)

    for alg in ['SAMME', 'SAMME.R']:
        clf = AdaBoostClassifier(algorithm=alg)

        clf.fit(X, y)
        importances = clf.feature_importances_

        assert_equal(importances.shape[0], 10)
        assert_equal((importances[:3, np.newaxis] >= importances[3:]).all(),
                     True)


def test_error():
    """Test that it gives proper exception on deficient input."""
    assert_raises(ValueError,
                  AdaBoostClassifier(learning_rate=-1).fit,
                  X, y_class)

    assert_raises(ValueError,
                  AdaBoostClassifier(algorithm="foo").fit,
                  X, y_class)

    assert_raises(TypeError,
                  AdaBoostClassifier(base_estimator=DummyRegressor()).fit,
                  X, y_class)

    assert_raises(TypeError,
                  AdaBoostRegressor(base_estimator=DummyClassifier(),
                                    random_state=0).fit,
                  X, y_regr)


def test_base_estimator():
    """Test different base estimators."""
    from sklearn.ensemble import RandomForestClassifier
    from sklearn.svm import SVC

    # XXX doesn't work with y_class because RF doesn't support classes_
    # Shouldn't AdaBoost run a LabelBinarizer?
    clf = AdaBoostClassifier(RandomForestClassifier())
    clf.fit(X, y_regr)

    clf = AdaBoostClassifier(SVC(), algorithm="SAMME")
    clf.fit(X, y_class)

    from sklearn.ensemble import RandomForestRegressor
    from sklearn.svm import SVR

    clf = AdaBoostRegressor(RandomForestRegressor(), random_state=0)
    clf.fit(X, y_regr)

    clf = AdaBoostRegressor(SVR(), random_state=0)
    clf.fit(X, y_regr)


if __name__ == "__main__":
    import nose
    nose.runmodule()
