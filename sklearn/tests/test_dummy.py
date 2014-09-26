from __future__ import division
import warnings
import numpy as np
import scipy.sparse as sp

from sklearn.base import clone
from sklearn.utils.testing import assert_array_equal
from sklearn.utils.testing import assert_array_almost_equal
from sklearn.utils.testing import assert_equal
from sklearn.utils.testing import assert_almost_equal
from sklearn.utils.testing import assert_raises
from sklearn.utils.testing import assert_true
from sklearn.utils.testing import assert_warns_message

from sklearn.dummy import DummyClassifier, DummyRegressor


def _check_predict_proba(clf, X, y):
    proba = clf.predict_proba(X)
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        # We know that we can have division by zero
        log_proba = clf.predict_log_proba(X)

    y = np.atleast_1d(y)
    if y.ndim == 1:
        y = np.reshape(y, (-1, 1))

    n_outputs = y.shape[1]
    n_samples = len(X)

    if n_outputs == 1:
        proba = [proba]
        log_proba = [log_proba]

    for k in range(n_outputs):
        assert_equal(proba[k].shape[0], n_samples)
        assert_equal(proba[k].shape[1], len(np.unique(y[:, k])))
        assert_array_equal(proba[k].sum(axis=1), np.ones(len(X)))
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            # We know that we can have division by zero
            assert_array_equal(np.log(proba[k]), log_proba[k])


def _check_behavior_2d(clf):
    # 1d case
    X = np.array([[0], [0], [0], [0]])  # ignored
    y = np.array([1, 2, 1, 1])
    est = clone(clf)
    est.fit(X, y)
    y_pred = est.predict(X)
    assert_equal(y.shape, y_pred.shape)

    # 2d case
    y = np.array([[1, 0],
                  [2, 0],
                  [1, 0],
                  [1, 3]])
    est = clone(clf)
    est.fit(X, y)
    y_pred = est.predict(X)
    assert_equal(y.shape, y_pred.shape)


def _check_behavior_2d_for_constant(clf):
    # 2d case only
    X = np.array([[0], [0], [0], [0]])  # ignored
    y = np.array([[1, 0, 5, 4, 3],
                  [2, 0, 1, 2, 5],
                  [1, 0, 4, 5, 2],
                  [1, 3, 3, 2, 0]])
    est = clone(clf)
    est.fit(X, y)
    y_pred = est.predict(X)
    assert_equal(y.shape, y_pred.shape)


def _check_equality_regressor(statistic, y_learn, y_pred_learn,
                              y_test, y_pred_test):
    assert_array_equal(np.tile(statistic, (y_learn.shape[0], 1)),
                       y_pred_learn)
    assert_array_equal(np.tile(statistic, (y_test.shape[0], 1)),
                       y_pred_test)


def test_most_frequent_strategy():
    X = [[0], [0], [0], [0]]  # ignored
    y = [1, 2, 1, 1]

    clf = DummyClassifier(strategy="most_frequent", random_state=0)
    clf.fit(X, y)
    assert_array_equal(clf.predict(X), np.ones(len(X)))
    _check_predict_proba(clf, X, y)


def test_most_frequent_strategy_multioutput():
    X = [[0], [0], [0], [0]]  # ignored
    y = np.array([[1, 0],
                  [2, 0],
                  [1, 0],
                  [1, 3]])

    n_samples = len(X)

    clf = DummyClassifier(strategy="most_frequent", random_state=0)
    clf.fit(X, y)
    assert_array_equal(clf.predict(X),
                       np.hstack([np.ones((n_samples, 1)),
                                  np.zeros((n_samples, 1))]))
    _check_predict_proba(clf, X, y)
    _check_behavior_2d(clf)


def test_stratified_strategy():
    X = [[0]] * 5  # ignored
    y = [1, 2, 1, 1, 2]
    clf = DummyClassifier(strategy="stratified", random_state=0)
    clf.fit(X, y)

    X = [[0]] * 500
    y_pred = clf.predict(X)
    p = np.bincount(y_pred) / float(len(X))
    assert_almost_equal(p[1], 3. / 5, decimal=1)
    assert_almost_equal(p[2], 2. / 5, decimal=1)
    _check_predict_proba(clf, X, y)


def test_stratified_strategy_multioutput():
    X = [[0]] * 5  # ignored
    y = np.array([[2, 1],
                  [2, 2],
                  [1, 1],
                  [1, 2],
                  [1, 1]])

    clf = DummyClassifier(strategy="stratified", random_state=0)
    clf.fit(X, y)

    X = [[0]] * 500
    y_pred = clf.predict(X)

    for k in range(y.shape[1]):
        p = np.bincount(y_pred[:, k]) / float(len(X))
        assert_almost_equal(p[1], 3. / 5, decimal=1)
        assert_almost_equal(p[2], 2. / 5, decimal=1)
        _check_predict_proba(clf, X, y)

    _check_behavior_2d(clf)


def test_uniform_strategy():
    X = [[0]] * 4  # ignored
    y = [1, 2, 1, 1]
    clf = DummyClassifier(strategy="uniform", random_state=0)
    clf.fit(X, y)

    X = [[0]] * 500
    y_pred = clf.predict(X)
    p = np.bincount(y_pred) / float(len(X))
    assert_almost_equal(p[1], 0.5, decimal=1)
    assert_almost_equal(p[2], 0.5, decimal=1)
    _check_predict_proba(clf, X, y)


def test_uniform_strategy_multioutput():
    X = [[0]] * 4  # ignored
    y = np.array([[2, 1],
                  [2, 2],
                  [1, 2],
                  [1, 1]])
    clf = DummyClassifier(strategy="uniform", random_state=0)
    clf.fit(X, y)

    X = [[0]] * 500
    y_pred = clf.predict(X)

    for k in range(y.shape[1]):
        p = np.bincount(y_pred[:, k]) / float(len(X))
        assert_almost_equal(p[1], 0.5, decimal=1)
        assert_almost_equal(p[2], 0.5, decimal=1)
        _check_predict_proba(clf, X, y)

    _check_behavior_2d(clf)


def test_string_labels():
    X = [[0]] * 5
    y = ["paris", "paris", "tokyo", "amsterdam", "berlin"]
    clf = DummyClassifier(strategy="most_frequent")
    clf.fit(X, y)
    assert_array_equal(clf.predict(X), ["paris"] * 5)


def test_classifier_exceptions():
    clf = DummyClassifier(strategy="unknown")
    assert_raises(ValueError, clf.fit, [], [])

    assert_raises(ValueError, clf.predict, [])
    assert_raises(ValueError, clf.predict_proba, [])


def test_mean_strategy_regressor():

    random_state = np.random.RandomState(seed=1)

    X = [[0]] * 4  # ignored
    y = random_state.randn(4)

    reg = DummyRegressor()
    reg.fit(X, y)
    assert_array_equal(reg.predict(X), [np.mean(y)] * len(X))


def test_mean_strategy_multioutput_regressor():

    random_state = np.random.RandomState(seed=1)

    X_learn = random_state.randn(10, 10)
    y_learn = random_state.randn(10, 5)

    mean = np.mean(y_learn, axis=0).reshape((1, -1))

    X_test = random_state.randn(20, 10)
    y_test = random_state.randn(20, 5)

    # Correctness oracle
    est = DummyRegressor()
    est.fit(X_learn, y_learn)
    y_pred_learn = est.predict(X_learn)
    y_pred_test = est.predict(X_test)

    _check_equality_regressor(mean, y_learn, y_pred_learn, y_test, y_pred_test)
    _check_behavior_2d(est)


def test_regressor_exceptions():
    reg = DummyRegressor()
    assert_raises(ValueError, reg.predict, [])


def test_median_strategy_regressor():

    random_state = np.random.RandomState(seed=1)

    X = [[0]] * 5  # ignored
    y = random_state.randn(5)

    reg = DummyRegressor(strategy="median")
    reg.fit(X, y)
    assert_array_equal(reg.predict(X), [np.median(y)] * len(X))


def test_median_strategy_multioutput_regressor():

    random_state = np.random.RandomState(seed=1)

    X_learn = random_state.randn(10, 10)
    y_learn = random_state.randn(10, 5)

    median = np.median(y_learn, axis=0).reshape((1, -1))

    X_test = random_state.randn(20, 10)
    y_test = random_state.randn(20, 5)

    # Correctness oracle
    est = DummyRegressor(strategy="median")
    est.fit(X_learn, y_learn)
    y_pred_learn = est.predict(X_learn)
    y_pred_test = est.predict(X_test)

    _check_equality_regressor(
        median, y_learn, y_pred_learn, y_test, y_pred_test)
    _check_behavior_2d(est)


def test_constant_strategy_regressor():

    random_state = np.random.RandomState(seed=1)

    X = [[0]] * 5  # ignored
    y = random_state.randn(5)

    reg = DummyRegressor(strategy="constant", constant=[43])
    reg.fit(X, y)
    assert_array_equal(reg.predict(X), [43] * len(X))

    reg = DummyRegressor(strategy="constant", constant=43)
    reg.fit(X, y)
    assert_array_equal(reg.predict(X), [43] * len(X))


def test_constant_strategy_multioutput_regressor():

    random_state = np.random.RandomState(seed=1)

    X_learn = random_state.randn(10, 10)
    y_learn = random_state.randn(10, 5)

    # test with 2d array
    constants = random_state.randn(5)

    X_test = random_state.randn(20, 10)
    y_test = random_state.randn(20, 5)

    # Correctness oracle
    est = DummyRegressor(strategy="constant", constant=constants)
    est.fit(X_learn, y_learn)
    y_pred_learn = est.predict(X_learn)
    y_pred_test = est.predict(X_test)

    _check_equality_regressor(
        constants, y_learn, y_pred_learn, y_test, y_pred_test)
    _check_behavior_2d_for_constant(est)


def test_y_mean_attribute_regressor():
    X = [[0]] * 5
    y = [1, 2, 4, 6, 8]
    # when strategy = 'mean'
    est = DummyRegressor(strategy='mean')
    est.fit(X, y)

    assert_equal(est.constant_, np.mean(y))


def test_unknown_strategey_regressor():
    X = [[0]] * 5
    y = [1, 2, 4, 6, 8]

    est = DummyRegressor(strategy='gona')
    assert_raises(ValueError, est.fit, X, y)


def test_constants_not_specified_regressor():
    X = [[0]] * 5
    y = [1, 2, 4, 6, 8]

    est = DummyRegressor(strategy='constant')
    assert_raises(TypeError, est.fit, X, y)


def test_constant_size_multioutput_regressor():
    random_state = np.random.RandomState(seed=1)
    X = random_state.randn(10, 10)
    y = random_state.randn(10, 5)

    est = DummyRegressor(strategy='constant', constant=[1, 2, 3, 4])
    assert_raises(ValueError, est.fit, X, y)


def test_constant_strategy():
    X = [[0], [0], [0], [0]]  # ignored
    y = [2, 1, 2, 2]

    clf = DummyClassifier(strategy="constant", random_state=0, constant=1)
    clf.fit(X, y)
    assert_array_equal(clf.predict(X), np.ones(len(X)))
    _check_predict_proba(clf, X, y)

    X = [[0], [0], [0], [0]]  # ignored
    y = ['two', 'one', 'two', 'two']
    clf = DummyClassifier(strategy="constant", random_state=0, constant='one')
    clf.fit(X, y)
    assert_array_equal(clf.predict(X), np.array(['one'] * 4))
    _check_predict_proba(clf, X, y)


def test_constant_strategy_multioutput():
    X = [[0], [0], [0], [0]]  # ignored
    y = np.array([[2, 3],
                  [1, 3],
                  [2, 3],
                  [2, 0]])

    n_samples = len(X)

    clf = DummyClassifier(strategy="constant", random_state=0,
                          constant=[1, 0])
    clf.fit(X, y)
    assert_array_equal(clf.predict(X),
                       np.hstack([np.ones((n_samples, 1)),
                                  np.zeros((n_samples, 1))]))
    _check_predict_proba(clf, X, y)


def test_constant_strategy_exceptions():
    X = [[0], [0], [0], [0]]  # ignored
    y = [2, 1, 2, 2]
    clf = DummyClassifier(strategy="constant", random_state=0)
    assert_raises(ValueError, clf.fit, X, y)
    clf = DummyClassifier(strategy="constant", random_state=0,
                          constant=[2, 0])
    assert_raises(ValueError, clf.fit, X, y)


def test_classification_sample_weight():
    X = [[0], [0], [1]]
    y = [0, 1, 0]
    sample_weight = [0.1, 1., 0.1]

    clf = DummyClassifier().fit(X, y, sample_weight)
    assert_array_almost_equal(clf.class_prior_, [0.2 / 1.2, 1. / 1.2])


def test_constant_strategy_sparse_target():
    X = [[0]] * 5  # ignored
    y = sp.csc_matrix(np.array([[0, 1],
                                [4, 0],
                                [1, 1],
                                [1, 4],
                                [1, 1]]))

    n_samples = len(X)

    clf = DummyClassifier(strategy="constant", random_state=0, constant=[1, 0])
    clf.fit(X, y)
    y_pred = clf.predict(X)
    assert_true(sp.issparse(y_pred))
    assert_array_equal(y_pred.toarray(), np.hstack([np.ones((n_samples, 1)),
                                                    np.zeros((n_samples, 1))]))


def test_uniform_strategy_sparse_target_warning():
    X = [[0]] * 5  # ignored
    y = sp.csc_matrix(np.array([[2, 1],
                                [2, 2],
                                [1, 4],
                                [4, 2],
                                [1, 1]]))

    clf = DummyClassifier(strategy="uniform", random_state=0)
    assert_warns_message(UserWarning,
                         "the uniform strategy would not save memory",
                         clf.fit, X, y)

    X = [[0]] * 500
    y_pred = clf.predict(X)

    for k in range(y.shape[1]):
        p = np.bincount(y_pred[:, k]) / float(len(X))
        assert_almost_equal(p[1], 1/3, decimal=1)
        assert_almost_equal(p[2], 1/3, decimal=1)
        assert_almost_equal(p[4], 1/3, decimal=1)


def test_stratified_strategy_sparse_target():
    X = [[0]] * 5  # ignored
    y = sp.csc_matrix(np.array([[4, 1],
                                [0, 0],
                                [1, 1],
                                [1, 4],
                                [1, 1]]))

    clf = DummyClassifier(strategy="stratified", random_state=0)
    clf.fit(X, y)

    X = [[0]] * 500
    y_pred = clf.predict(X)
    assert_true(sp.issparse(y_pred))
    y_pred = y_pred.toarray()

    for k in range(y.shape[1]):
        p = np.bincount(y_pred[:, k]) / float(len(X))
        assert_almost_equal(p[1], 3. / 5, decimal=1)
        assert_almost_equal(p[0], 1. / 5, decimal=1)
        assert_almost_equal(p[4], 1. / 5, decimal=1)


def test_most_frequent_strategy_sparse_target():
    X = [[0]] * 5  # ignored
    y = sp.csc_matrix(np.array([[1, 0],
                                [1, 3],
                                [4, 0],
                                [0, 1],
                                [1, 0]]))

    n_samples = len(X)
    clf = DummyClassifier(strategy="most_frequent", random_state=0)
    clf.fit(X, y)

    y_pred = clf.predict(X)
    assert_true(sp.issparse(y_pred))
    assert_array_equal(y_pred.toarray(), np.hstack([np.ones((n_samples, 1)),
                                                    np.zeros((n_samples, 1))]))


if __name__ == '__main__':
    import nose
    nose.runmodule()
