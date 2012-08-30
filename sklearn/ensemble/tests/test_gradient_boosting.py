"""
Testing for the gradient boosting module (sklearn.ensemble.gradient_boosting).
"""

import numpy as np
from numpy.testing import assert_array_equal
from numpy.testing import assert_array_almost_equal
from numpy.testing import assert_equal

from nose.tools import assert_raises

from sklearn.metrics import mean_squared_error
from sklearn.utils import check_random_state

from sklearn.ensemble import GradientBoostingClassifier
from sklearn.ensemble import GradientBoostingRegressor

from sklearn import datasets

# toy sample
X = [[-2, -1], [-1, -1], [-1, -2], [1, 1], [1, 2], [2, 1]]
y = [-1, -1, -1, 1, 1, 1]
T = [[-1, -1], [2, 2], [3, 2]]
true_result = [-1, 1, 1]

rng = np.random.RandomState(0)
# also load the boston dataset
# and randomly permute it
boston = datasets.load_boston()
perm = rng.permutation(boston.target.size)
boston.data = boston.data[perm]
boston.target = boston.target[perm]

# also load the iris dataset
# and randomly permute it
iris = datasets.load_iris()
perm = rng.permutation(iris.target.size)
iris.data = iris.data[perm]
iris.target = iris.target[perm]


def test_classification_toy():
    """Check classification on a toy dataset."""
    clf = GradientBoostingClassifier(n_estimators=100, random_state=1)

    assert_raises(ValueError, clf.predict, T)

    clf.fit(X, y)
    assert_array_equal(clf.predict(T), true_result)
    assert_equal(100, len(clf.estimators_))

    deviance_decrease = (clf.train_score_[:-1] - clf.train_score_[1:])
    assert np.any(deviance_decrease >= 0.0), \
           "Train deviance does not monotonically decrease."


def test_parameter_checks():
    """Check input parameter validation."""
    assert_raises(ValueError, GradientBoostingClassifier, n_estimators=0)
    assert_raises(ValueError, GradientBoostingClassifier, n_estimators=-1)

    assert_raises(ValueError, GradientBoostingClassifier, learn_rate=0.0)
    assert_raises(ValueError, GradientBoostingClassifier, learn_rate=-1.0)

    assert_raises(ValueError, GradientBoostingRegressor, loss='foobar')

    assert_raises(ValueError, GradientBoostingClassifier,
                  min_samples_split=0.0)
    assert_raises(ValueError, GradientBoostingClassifier,
                  min_samples_split=-1.0)

    assert_raises(ValueError, GradientBoostingClassifier, min_samples_leaf=0)
    assert_raises(ValueError, GradientBoostingClassifier, min_samples_leaf=-1.)

    assert_raises(ValueError, GradientBoostingClassifier, subsample=0.0)
    assert_raises(ValueError, GradientBoostingClassifier, subsample=1.1)
    assert_raises(ValueError, GradientBoostingClassifier, subsample=-0.1)

    assert_raises(ValueError, GradientBoostingClassifier, max_depth=-0.1)
    assert_raises(ValueError, GradientBoostingClassifier, max_depth=0)

    assert_raises(ValueError, GradientBoostingClassifier, init={})

    # test fit before feature importance
    assert_raises(ValueError,
                  lambda: GradientBoostingClassifier().feature_importances_)

    # binomial deviance requires ``n_classes == 2``.
    assert_raises(ValueError,
                  lambda X, y: GradientBoostingClassifier(
                      loss='bdeviance').fit(X, y),
                  X, [0, 0, 1, 1, 2, 2])

    # multinomial deviance requires ``n_classes > 2``.
    assert_raises(ValueError,
                  lambda X, y: GradientBoostingClassifier(
                      loss='mdeviance').fit(X, y),
                  X, [0, 0, 1, 1, 1, 0])

    # deviance requires ``n_classes >= 2``.
    assert_raises(ValueError,
                  lambda X, y: GradientBoostingClassifier(
                      loss='deviance').fit(X, y),
                  X, [0, 0, 0, 0])


def test_classification_synthetic():
    """Test GradientBoostingClassifier on synthetic dataset used by
    Hastie et al. in ESLII Example 12.7. """
    X, y = datasets.make_hastie_10_2(n_samples=12000, random_state=1)

    X_train, X_test = X[:2000], X[2000:]
    y_train, y_test = y[:2000], y[2000:]

    gbrt = GradientBoostingClassifier(n_estimators=100, min_samples_split=1,
                                      max_depth=1,
                                      learn_rate=1.0, random_state=0)
    gbrt.fit(X_train, y_train)
    error_rate = (1.0 - gbrt.score(X_test, y_test))
    assert error_rate < 0.085, \
           "GB failed with error %.4f" % error_rate

    gbrt = GradientBoostingClassifier(n_estimators=200, min_samples_split=1,
                                      max_depth=1,
                                      learn_rate=1.0, subsample=0.5,
                                      random_state=0)
    gbrt.fit(X_train, y_train)
    error_rate = (1.0 - gbrt.score(X_test, y_test))
    assert error_rate < 0.08, \
           "Stochastic GB failed with error %.4f" % error_rate


def test_boston():
    """Check consistency on dataset boston house prices with least squares
    and least absolute deviation. """
    for loss in ("ls", "lad", "huber"):
        clf = GradientBoostingRegressor(n_estimators=100, loss=loss,
                                        max_depth=4,
                                        min_samples_split=1, random_state=1)
        assert_raises(ValueError, clf.predict, boston.data)
        clf.fit(boston.data, boston.target)
        y_pred = clf.predict(boston.data)
        mse = mean_squared_error(boston.target, y_pred)
        assert mse < 6.0, "Failed with loss %s and mse = %.4f" % (loss, mse)


def test_iris():
    """Check consistency on dataset iris."""
    for subsample in (1.0, 0.5):
        clf = GradientBoostingClassifier(n_estimators=100, loss='deviance',
                                         random_state=1, subsample=subsample)
        clf.fit(iris.data, iris.target)
        score = clf.score(iris.data, iris.target)
        assert score > 0.9, "Failed with subsample %.1f " \
               "and score = %f" % (subsample, score)


def test_regression_synthetic():
    """Test on synthetic regression datasets used in Leo Breiman,
    `Bagging Predictors?. Machine Learning 24(2): 123-140 (1996). """
    random_state = check_random_state(1)
    regression_params = {'n_estimators': 100, 'max_depth': 4,
                         'min_samples_split': 1, 'learn_rate': 0.1,
                         'loss': 'ls'}

    # Friedman1
    X, y = datasets.make_friedman1(n_samples=1200,
                                   random_state=random_state, noise=1.0)
    X_train, y_train = X[:200], y[:200]
    X_test, y_test = X[200:], y[200:]
    clf = GradientBoostingRegressor()
    clf.fit(X_train, y_train)
    mse = mean_squared_error(y_test, clf.predict(X_test))
    assert mse < 5.0, "Failed on Friedman1 with mse = %.4f" % mse

    # Friedman2
    X, y = datasets.make_friedman2(n_samples=1200, random_state=random_state)
    X_train, y_train = X[:200], y[:200]
    X_test, y_test = X[200:], y[200:]
    clf = GradientBoostingRegressor(**regression_params)
    clf.fit(X_train, y_train)
    mse = mean_squared_error(y_test, clf.predict(X_test))
    assert mse < 1700.0, "Failed on Friedman2 with mse = %.4f" % mse

    # Friedman3
    X, y = datasets.make_friedman3(n_samples=1200, random_state=random_state)
    X_train, y_train = X[:200], y[:200]
    X_test, y_test = X[200:], y[200:]
    clf = GradientBoostingRegressor(**regression_params)
    clf.fit(X_train, y_train)
    mse = mean_squared_error(y_test, clf.predict(X_test))
    assert mse < 0.015, "Failed on Friedman3 with mse = %.4f" % mse


# def test_feature_importances():
#     X = np.array(boston.data, dtype=np.float32)
#     y = np.array(boston.target, dtype=np.float32)

#     clf = GradientBoostingRegressor(n_estimators=100, max_depth=5,
#                                     min_samples_split=1, random_state=1)
#     clf.fit(X, y)
#     feature_importances = clf.feature_importances_

#     # true feature importance ranking
#     true_ranking = np.array([3, 1, 8, 2, 10, 9, 4, 11, 0, 6, 7, 5, 12])

#     assert_array_equal(true_ranking, feature_importances.argsort())


def test_probability():
    """Predict probabilities."""
    clf = GradientBoostingClassifier(n_estimators=100, random_state=1)

    assert_raises(ValueError, clf.predict_proba, T)

    clf.fit(X, y)
    assert_array_equal(clf.predict(T), true_result)

    # check if probabilities are in [0, 1].
    y_proba = clf.predict_proba(T)
    assert np.all(y_proba >= 0.0)
    assert np.all(y_proba <= 1.0)

    # derive predictions from probabilities
    y_pred = clf.classes_.take(y_proba.argmax(axis=1), axis=0)
    assert_array_equal(y_pred, true_result)


def test_check_inputs():
    """Test input checks (shape and type of X and y)."""
    clf = GradientBoostingClassifier(n_estimators=100, random_state=1)
    assert_raises(ValueError, clf.fit, X, y + [0, 1])

    from scipy import sparse
    X_sparse = sparse.csr_matrix(X)
    clf = GradientBoostingClassifier(n_estimators=100, random_state=1)
    assert_raises(TypeError, clf.fit, X_sparse, y)

    clf = GradientBoostingClassifier().fit(X, y)
    assert_raises(TypeError, clf.predict, X_sparse)


def test_check_inputs_predict():
    """X has wrong shape """
    clf = GradientBoostingClassifier(n_estimators=100, random_state=1)
    clf.fit(X, y)

    x = np.array([1.0, 2.0])[:, np.newaxis]
    assert_raises(ValueError, clf.predict, x)

    x = np.array([])
    assert_raises(ValueError, clf.predict, x)

    x = np.array([1.0, 2.0, 3.0])[:, np.newaxis]
    assert_raises(ValueError, clf.predict, x)

    clf = GradientBoostingRegressor(n_estimators=100, random_state=1)
    clf.fit(X, rng.rand(len(X)))

    x = np.array([1.0, 2.0])[:, np.newaxis]
    assert_raises(ValueError, clf.predict, x)

    x = np.array([])
    assert_raises(ValueError, clf.predict, x)

    x = np.array([1.0, 2.0, 3.0])[:, np.newaxis]
    assert_raises(ValueError, clf.predict, x)


def test_check_max_features():
    """test if max_features is valid. """
    clf = GradientBoostingRegressor(n_estimators=100, random_state=1,
                                    max_features=0)
    assert_raises(ValueError, clf.fit, X, y)

    clf = GradientBoostingRegressor(n_estimators=100, random_state=1,
                                    max_features=(len(X[0]) + 1))
    assert_raises(ValueError, clf.fit, X, y)


def test_staged_predict():
    """Test whether staged decision function eventually gives
    the same prediction.
    """
    X, y = datasets.make_friedman1(n_samples=1200,
                                   random_state=1, noise=1.0)
    X_train, y_train = X[:200], y[:200]
    X_test, y_test = X[200:], y[200:]
    clf = GradientBoostingRegressor()
    # test raise ValueError if not fitted
    assert_raises(ValueError, lambda X: np.fromiter(
        clf.staged_predict(X), dtype=np.float64), X_test)

    clf.fit(X_train, y_train)
    y_pred = clf.predict(X_test)

    # test if prediction for last stage equals ``predict``
    for y in clf.staged_predict(X_test):
        assert_equal(y.shape, y_pred.shape)

    assert_array_equal(y_pred, y)


def test_serialization():
    """Check model serialization."""
    clf = GradientBoostingClassifier(n_estimators=100, random_state=1)

    clf.fit(X, y)
    assert_array_equal(clf.predict(T), true_result)
    assert_equal(100, len(clf.estimators_))

    try:
        import cPickle as pickle
    except ImportError:
        import pickle

    serialized_clf = pickle.dumps(clf, protocol=pickle.HIGHEST_PROTOCOL)
    clf = None
    clf = pickle.loads(serialized_clf)
    assert_array_equal(clf.predict(T), true_result)
    assert_equal(100, len(clf.estimators_))


def test_degenerate_targets():
    """Check if we can fit even though all targets are equal. """
    clf = GradientBoostingClassifier(n_estimators=100, random_state=1)

    # classifier should raise exception
    assert_raises(ValueError, clf.fit, X, np.ones(len(X)))

    clf = GradientBoostingRegressor(n_estimators=100, random_state=1)
    clf.fit(X, np.ones(len(X)))
    clf.predict(rng.rand(2))
    assert_array_equal(np.ones((1,), dtype=np.float64),
                       clf.predict(rng.rand(2)))


def test_quantile_loss():
    """Check if quantile loss with alpha=0.5 equals lad. """
    clf_quantile = GradientBoostingRegressor(n_estimators=100, loss='quantile',
                                             max_depth=4, alpha=0.5,
                                             random_state=7)

    clf_quantile.fit(boston.data, boston.target)
    y_quantile = clf_quantile.predict(boston.data)

    clf_lad = GradientBoostingRegressor(n_estimators=100, loss='lad',
                                        max_depth=4, random_state=7)

    clf_lad.fit(boston.data, boston.target)
    y_lad = clf_lad.predict(boston.data)
    assert_array_almost_equal(y_quantile, y_lad, decimal=4)


def test_symbol_labels():
    """Test with non-integer class labels. """
    clf = GradientBoostingClassifier(n_estimators=100, random_state=1)

    symbol_y = map(str, y)

    clf.fit(X, symbol_y)
    assert_array_equal(clf.predict(T), map(str, true_result))
    assert_equal(100, len(clf.estimators_))


def test_float_class_labels():
    """Test with float class labels. """
    clf = GradientBoostingClassifier(n_estimators=100, random_state=1)

    float_y = np.asarray(y, dtype=np.float32)

    clf.fit(X, float_y)
    assert_array_equal(clf.predict(T),
                       np.asarray(true_result, dtype=np.float32))
    assert_equal(100, len(clf.estimators_))


def test_shape_y():
    """Test with float class labels. """
    clf = GradientBoostingClassifier(n_estimators=100, random_state=1)

    y_ = np.asarray(y, dtype=np.int32)
    y_ = y_[:, np.newaxis]

    clf.fit(X, y_)
    assert_array_equal(clf.predict(T), true_result)
    assert_equal(100, len(clf.estimators_))


def test_mem_layout():
    """Test with different memory layouts of X and y"""
    X_ = np.asfortranarray(X)
    clf = GradientBoostingClassifier(n_estimators=100, random_state=1)
    clf.fit(X_, y)
    assert_array_equal(clf.predict(T), true_result)
    assert_equal(100, len(clf.estimators_))

    X_ = np.ascontiguousarray(X)
    clf = GradientBoostingClassifier(n_estimators=100, random_state=1)
    clf.fit(X_, y)
    assert_array_equal(clf.predict(T), true_result)
    assert_equal(100, len(clf.estimators_))

    y_ = np.asarray(y, dtype=np.int32)
    y_ = y_[:, np.newaxis]
    y_ = np.ascontiguousarray(y_)
    clf = GradientBoostingClassifier(n_estimators=100, random_state=1)
    clf.fit(X, y_)
    assert_array_equal(clf.predict(T), true_result)
    assert_equal(100, len(clf.estimators_))

    y_ = np.asarray(y, dtype=np.int32)
    y_ = y_[:, np.newaxis]
    y_ = np.asfortranarray(y_)
    clf = GradientBoostingClassifier(n_estimators=100, random_state=1)
    clf.fit(X, y_)
    assert_array_equal(clf.predict(T), true_result)
    assert_equal(100, len(clf.estimators_))
