"""
Testing for the gradient boosting module (sklearn.ensemble.gradient_boosting).
"""
import warnings
import numpy as np

from itertools import product

from scipy.sparse import csr_matrix
from scipy.sparse import csc_matrix
from scipy.sparse import coo_matrix

from sklearn import datasets
from sklearn.base import clone
from sklearn.ensemble import GradientBoostingClassifier
from sklearn.ensemble import GradientBoostingRegressor
from sklearn.ensemble.gradient_boosting import ZeroEstimator
from sklearn.metrics import mean_squared_error
from sklearn.utils import check_random_state, tosequence
from sklearn.utils.testing import assert_almost_equal
from sklearn.utils.testing import assert_array_almost_equal
from sklearn.utils.testing import assert_array_equal
from sklearn.utils.testing import assert_equal
from sklearn.utils.testing import assert_greater
from sklearn.utils.testing import assert_less
from sklearn.utils.testing import assert_raises
from sklearn.utils.testing import assert_true
from sklearn.utils.testing import assert_warns
from sklearn.utils.testing import assert_warns_message
from sklearn.utils.testing import skip_if_32bit
from sklearn.exceptions import DataConversionWarning
from sklearn.exceptions import NotFittedError

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


def check_classification_toy(presort, loss):
    # Check classification on a toy dataset.
    clf = GradientBoostingClassifier(loss=loss, n_estimators=10,
                                     random_state=1, presort=presort)

    assert_raises(ValueError, clf.predict, T)

    clf.fit(X, y)
    assert_array_equal(clf.predict(T), true_result)
    assert_equal(10, len(clf.estimators_))

    deviance_decrease = (clf.train_score_[:-1] - clf.train_score_[1:])
    assert_true(np.any(deviance_decrease >= 0.0))

    leaves = clf.apply(X)
    assert_equal(leaves.shape, (6, 10, 1))


def test_classification_toy():
    for presort, loss in product(('auto', True, False),
                                 ('deviance', 'exponential')):
        yield check_classification_toy, presort, loss


def test_parameter_checks():
    # Check input parameter validation.

    assert_raises(ValueError,
                  GradientBoostingClassifier(n_estimators=0).fit, X, y)
    assert_raises(ValueError,
                  GradientBoostingClassifier(n_estimators=-1).fit, X, y)

    assert_raises(ValueError,
                  GradientBoostingClassifier(learning_rate=0.0).fit, X, y)
    assert_raises(ValueError,
                  GradientBoostingClassifier(learning_rate=-1.0).fit, X, y)

    assert_raises(ValueError,
                  GradientBoostingClassifier(loss='foobar').fit, X, y)

    assert_raises(ValueError,
                  GradientBoostingClassifier(min_samples_split=0.0).fit, X, y)
    assert_raises(ValueError,
                  GradientBoostingClassifier(min_samples_split=-1.0).fit, X, y)
    assert_raises(ValueError,
                  GradientBoostingClassifier(min_samples_split=1.1).fit, X, y)

    assert_raises(ValueError,
                  GradientBoostingClassifier(min_samples_leaf=0).fit, X, y)
    assert_raises(ValueError,
                  GradientBoostingClassifier(min_samples_leaf=-1.0).fit, X, y)

    assert_raises(ValueError,
                  GradientBoostingClassifier(min_weight_fraction_leaf=-1.).fit,
                  X, y)
    assert_raises(ValueError,
                  GradientBoostingClassifier(min_weight_fraction_leaf=0.6).fit,
                  X, y)

    assert_raises(ValueError,
                  GradientBoostingClassifier(subsample=0.0).fit, X, y)
    assert_raises(ValueError,
                  GradientBoostingClassifier(subsample=1.1).fit, X, y)
    assert_raises(ValueError,
                  GradientBoostingClassifier(subsample=-0.1).fit, X, y)

    assert_raises(ValueError,
                  GradientBoostingClassifier(max_depth=-0.1).fit, X, y)
    assert_raises(ValueError,
                  GradientBoostingClassifier(max_depth=0).fit, X, y)

    assert_raises(ValueError,
                  GradientBoostingClassifier(init={}).fit, X, y)

    # test fit before feature importance
    assert_raises(ValueError,
                  lambda: GradientBoostingClassifier().feature_importances_)

    # deviance requires ``n_classes >= 2``.
    assert_raises(ValueError,
                  lambda X, y: GradientBoostingClassifier(
                      loss='deviance').fit(X, y),
                  X, [0, 0, 0, 0])


def test_loss_function():
    assert_raises(ValueError,
                  GradientBoostingClassifier(loss='ls').fit, X, y)
    assert_raises(ValueError,
                  GradientBoostingClassifier(loss='lad').fit, X, y)
    assert_raises(ValueError,
                  GradientBoostingClassifier(loss='quantile').fit, X, y)
    assert_raises(ValueError,
                  GradientBoostingClassifier(loss='huber').fit, X, y)
    assert_raises(ValueError,
                  GradientBoostingRegressor(loss='deviance').fit, X, y)
    assert_raises(ValueError,
                  GradientBoostingRegressor(loss='exponential').fit, X, y)


def check_classification_synthetic(presort, loss):
    # Test GradientBoostingClassifier on synthetic dataset used by
    # Hastie et al. in ESLII Example 12.7.
    X, y = datasets.make_hastie_10_2(n_samples=12000, random_state=1)

    X_train, X_test = X[:2000], X[2000:]
    y_train, y_test = y[:2000], y[2000:]

    gbrt = GradientBoostingClassifier(n_estimators=100, min_samples_split=2,
                                      max_depth=1, loss=loss,
                                      learning_rate=1.0, random_state=0)
    gbrt.fit(X_train, y_train)
    error_rate = (1.0 - gbrt.score(X_test, y_test))
    assert_less(error_rate, 0.09)

    gbrt = GradientBoostingClassifier(n_estimators=200, min_samples_split=2,
                                      max_depth=1, loss=loss,
                                      learning_rate=1.0, subsample=0.5,
                                      random_state=0,
                                      presort=presort)
    gbrt.fit(X_train, y_train)
    error_rate = (1.0 - gbrt.score(X_test, y_test))
    assert_less(error_rate, 0.08)


def test_classification_synthetic():
    for presort, loss in product(('auto', True, False), ('deviance', 'exponential')):
        yield check_classification_synthetic, presort, loss


def check_boston(presort, loss, subsample):
    # Check consistency on dataset boston house prices with least squares
    # and least absolute deviation.
    ones = np.ones(len(boston.target))
    last_y_pred = None
    for sample_weight in None, ones, 2 * ones:
        clf = GradientBoostingRegressor(n_estimators=100,
                                        loss=loss,
                                        max_depth=4,
                                        subsample=subsample,
                                        min_samples_split=2,
                                        random_state=1,
                                        presort=presort)

        assert_raises(ValueError, clf.predict, boston.data)
        clf.fit(boston.data, boston.target,
                sample_weight=sample_weight)
        leaves = clf.apply(boston.data)
        assert_equal(leaves.shape, (506, 100))

        y_pred = clf.predict(boston.data)
        mse = mean_squared_error(boston.target, y_pred)
        assert_less(mse, 6.0)

        if last_y_pred is not None:
            assert_array_almost_equal(last_y_pred, y_pred)

        last_y_pred = y_pred


def test_boston():
    for presort, loss, subsample in product(('auto', True, False),
                                            ('ls', 'lad', 'huber'),
                                            (1.0, 0.5)):
        yield check_boston, presort, loss, subsample


def check_iris(presort, subsample, sample_weight):
    # Check consistency on dataset iris.
    clf = GradientBoostingClassifier(n_estimators=100,
                                     loss='deviance',
                                     random_state=1,
                                     subsample=subsample,
                                     presort=presort)
    clf.fit(iris.data, iris.target, sample_weight=sample_weight)
    score = clf.score(iris.data, iris.target)
    assert_greater(score, 0.9)

    leaves = clf.apply(iris.data)
    assert_equal(leaves.shape, (150, 100, 3))


def test_iris():
    ones = np.ones(len(iris.target))
    for presort, subsample, sample_weight in product(('auto', True, False),
                                                     (1.0, 0.5),
                                                     (None, ones)):
        yield check_iris, presort, subsample, sample_weight


def test_regression_synthetic():
    # Test on synthetic regression datasets used in Leo Breiman,
    # `Bagging Predictors?. Machine Learning 24(2): 123-140 (1996).
    random_state = check_random_state(1)
    regression_params = {'n_estimators': 100, 'max_depth': 4,
                         'min_samples_split': 2, 'learning_rate': 0.1,
                         'loss': 'ls'}

    # Friedman1
    X, y = datasets.make_friedman1(n_samples=1200,
                                   random_state=random_state,
                                   noise=1.0)
    X_train, y_train = X[:200], y[:200]
    X_test, y_test = X[200:], y[200:]

    for presort in True, False:
        clf = GradientBoostingRegressor(presort=presort)
        clf.fit(X_train, y_train)
        mse = mean_squared_error(y_test, clf.predict(X_test))
        assert_less(mse, 5.0)

    # Friedman2
    X, y = datasets.make_friedman2(n_samples=1200, random_state=random_state)
    X_train, y_train = X[:200], y[:200]
    X_test, y_test = X[200:], y[200:]

    for presort in True, False:
        regression_params['presort'] = presort
        clf = GradientBoostingRegressor(**regression_params)
        clf.fit(X_train, y_train)
        mse = mean_squared_error(y_test, clf.predict(X_test))
        assert_less(mse, 1700.0)

    # Friedman3
    X, y = datasets.make_friedman3(n_samples=1200, random_state=random_state)
    X_train, y_train = X[:200], y[:200]
    X_test, y_test = X[200:], y[200:]

    for presort in True, False:
        regression_params['presort'] = presort
        clf = GradientBoostingRegressor(**regression_params)
        clf.fit(X_train, y_train)
        mse = mean_squared_error(y_test, clf.predict(X_test))
        assert_less(mse, 0.015)


def test_feature_importances():
    X = np.array(boston.data, dtype=np.float32)
    y = np.array(boston.target, dtype=np.float32)

    for presort in True, False:
        clf = GradientBoostingRegressor(n_estimators=100, max_depth=5,
                                        min_samples_split=2, random_state=1,
                                        presort=presort)
        clf.fit(X, y)
        assert_true(hasattr(clf, 'feature_importances_'))


def test_probability_log():
    # Predict probabilities.
    clf = GradientBoostingClassifier(n_estimators=100, random_state=1)

    assert_raises(ValueError, clf.predict_proba, T)

    clf.fit(X, y)
    assert_array_equal(clf.predict(T), true_result)

    # check if probabilities are in [0, 1].
    y_proba = clf.predict_proba(T)
    assert_true(np.all(y_proba >= 0.0))
    assert_true(np.all(y_proba <= 1.0))

    # derive predictions from probabilities
    y_pred = clf.classes_.take(y_proba.argmax(axis=1), axis=0)
    assert_array_equal(y_pred, true_result)


def test_check_inputs():
    # Test input checks (shape and type of X and y).
    clf = GradientBoostingClassifier(n_estimators=100, random_state=1)
    assert_raises(ValueError, clf.fit, X, y + [0, 1])

    clf = GradientBoostingClassifier(n_estimators=100, random_state=1)
    assert_raises(ValueError, clf.fit, X, y,
                  sample_weight=([1] * len(y)) + [0, 1])


def test_check_inputs_predict():
    # X has wrong shape
    clf = GradientBoostingClassifier(n_estimators=100, random_state=1)
    clf.fit(X, y)

    x = np.array([1.0, 2.0])[:, np.newaxis]
    assert_raises(ValueError, clf.predict, x)

    x = np.array([[]])
    assert_raises(ValueError, clf.predict, x)

    x = np.array([1.0, 2.0, 3.0])[:, np.newaxis]
    assert_raises(ValueError, clf.predict, x)

    clf = GradientBoostingRegressor(n_estimators=100, random_state=1)
    clf.fit(X, rng.rand(len(X)))

    x = np.array([1.0, 2.0])[:, np.newaxis]
    assert_raises(ValueError, clf.predict, x)

    x = np.array([[]])
    assert_raises(ValueError, clf.predict, x)

    x = np.array([1.0, 2.0, 3.0])[:, np.newaxis]
    assert_raises(ValueError, clf.predict, x)


def test_check_max_features():
    # test if max_features is valid.
    clf = GradientBoostingRegressor(n_estimators=100, random_state=1,
                                    max_features=0)
    assert_raises(ValueError, clf.fit, X, y)

    clf = GradientBoostingRegressor(n_estimators=100, random_state=1,
                                    max_features=(len(X[0]) + 1))
    assert_raises(ValueError, clf.fit, X, y)

    clf = GradientBoostingRegressor(n_estimators=100, random_state=1,
                                    max_features=-0.1)
    assert_raises(ValueError, clf.fit, X, y)


def test_max_feature_regression():
    # Test to make sure random state is set properly.
    X, y = datasets.make_hastie_10_2(n_samples=12000, random_state=1)

    X_train, X_test = X[:2000], X[2000:]
    y_train, y_test = y[:2000], y[2000:]

    gbrt = GradientBoostingClassifier(n_estimators=100, min_samples_split=5,
                                      max_depth=2, learning_rate=.1,
                                      max_features=2, random_state=1)
    gbrt.fit(X_train, y_train)
    deviance = gbrt.loss_(y_test, gbrt.decision_function(X_test))
    assert_true(deviance < 0.5, "GB failed with deviance %.4f" % deviance)


def test_max_feature_auto():
    # Test if max features is set properly for floats and str.
    X, y = datasets.make_hastie_10_2(n_samples=12000, random_state=1)
    _, n_features = X.shape

    X_train = X[:2000]
    y_train = y[:2000]

    gbrt = GradientBoostingClassifier(n_estimators=1, max_features='auto')
    gbrt.fit(X_train, y_train)
    assert_equal(gbrt.max_features_, int(np.sqrt(n_features)))

    gbrt = GradientBoostingRegressor(n_estimators=1, max_features='auto')
    gbrt.fit(X_train, y_train)
    assert_equal(gbrt.max_features_, n_features)

    gbrt = GradientBoostingRegressor(n_estimators=1, max_features=0.3)
    gbrt.fit(X_train, y_train)
    assert_equal(gbrt.max_features_, int(n_features * 0.3))

    gbrt = GradientBoostingRegressor(n_estimators=1, max_features='sqrt')
    gbrt.fit(X_train, y_train)
    assert_equal(gbrt.max_features_, int(np.sqrt(n_features)))

    gbrt = GradientBoostingRegressor(n_estimators=1, max_features='log2')
    gbrt.fit(X_train, y_train)
    assert_equal(gbrt.max_features_, int(np.log2(n_features)))

    gbrt = GradientBoostingRegressor(n_estimators=1,
                                     max_features=0.01 / X.shape[1])
    gbrt.fit(X_train, y_train)
    assert_equal(gbrt.max_features_, 1)


def test_staged_predict():
    # Test whether staged decision function eventually gives
    # the same prediction.
    X, y = datasets.make_friedman1(n_samples=1200,
                                   random_state=1, noise=1.0)
    X_train, y_train = X[:200], y[:200]
    X_test = X[200:]
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


def test_staged_predict_proba():
    # Test whether staged predict proba eventually gives
    # the same prediction.
    X, y = datasets.make_hastie_10_2(n_samples=1200,
                                     random_state=1)
    X_train, y_train = X[:200], y[:200]
    X_test, y_test = X[200:], y[200:]
    clf = GradientBoostingClassifier(n_estimators=20)
    # test raise NotFittedError if not fitted
    assert_raises(NotFittedError, lambda X: np.fromiter(
        clf.staged_predict_proba(X), dtype=np.float64), X_test)

    clf.fit(X_train, y_train)

    # test if prediction for last stage equals ``predict``
    for y_pred in clf.staged_predict(X_test):
        assert_equal(y_test.shape, y_pred.shape)

    assert_array_equal(clf.predict(X_test), y_pred)

    # test if prediction for last stage equals ``predict_proba``
    for staged_proba in clf.staged_predict_proba(X_test):
        assert_equal(y_test.shape[0], staged_proba.shape[0])
        assert_equal(2, staged_proba.shape[1])

    assert_array_equal(clf.predict_proba(X_test), staged_proba)


def test_staged_functions_defensive():
    # test that staged_functions make defensive copies
    rng = np.random.RandomState(0)
    X = rng.uniform(size=(10, 3))
    y = (4 * X[:, 0]).astype(np.int) + 1  # don't predict zeros
    for estimator in [GradientBoostingRegressor(),
                      GradientBoostingClassifier()]:
        estimator.fit(X, y)
        for func in ['predict', 'decision_function', 'predict_proba']:
            staged_func = getattr(estimator, "staged_" + func, None)
            if staged_func is None:
                # regressor has no staged_predict_proba
                continue
            with warnings.catch_warnings(record=True):
                staged_result = list(staged_func(X))
            staged_result[1][:] = 0
            assert_true(np.all(staged_result[0] != 0))


def test_serialization():
    # Check model serialization.
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
    # Check if we can fit even though all targets are equal.
    clf = GradientBoostingClassifier(n_estimators=100, random_state=1)

    # classifier should raise exception
    assert_raises(ValueError, clf.fit, X, np.ones(len(X)))

    clf = GradientBoostingRegressor(n_estimators=100, random_state=1)
    clf.fit(X, np.ones(len(X)))
    clf.predict([rng.rand(2)])
    assert_array_equal(np.ones((1,), dtype=np.float64),
                       clf.predict([rng.rand(2)]))


def test_quantile_loss():
    # Check if quantile loss with alpha=0.5 equals lad.
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
    # Test with non-integer class labels.
    clf = GradientBoostingClassifier(n_estimators=100, random_state=1)

    symbol_y = tosequence(map(str, y))

    clf.fit(X, symbol_y)
    assert_array_equal(clf.predict(T), tosequence(map(str, true_result)))
    assert_equal(100, len(clf.estimators_))


def test_float_class_labels():
    # Test with float class labels.
    clf = GradientBoostingClassifier(n_estimators=100, random_state=1)

    float_y = np.asarray(y, dtype=np.float32)

    clf.fit(X, float_y)
    assert_array_equal(clf.predict(T),
                       np.asarray(true_result, dtype=np.float32))
    assert_equal(100, len(clf.estimators_))


def test_shape_y():
    # Test with float class labels.
    clf = GradientBoostingClassifier(n_estimators=100, random_state=1)

    y_ = np.asarray(y, dtype=np.int32)
    y_ = y_[:, np.newaxis]

    # This will raise a DataConversionWarning that we want to
    # "always" raise, elsewhere the warnings gets ignored in the
    # later tests, and the tests that check for this warning fail
    assert_warns(DataConversionWarning, clf.fit, X, y_)
    assert_array_equal(clf.predict(T), true_result)
    assert_equal(100, len(clf.estimators_))


def test_mem_layout():
    # Test with different memory layouts of X and y
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
    y_ = np.ascontiguousarray(y_)
    clf = GradientBoostingClassifier(n_estimators=100, random_state=1)
    clf.fit(X, y_)
    assert_array_equal(clf.predict(T), true_result)
    assert_equal(100, len(clf.estimators_))

    y_ = np.asarray(y, dtype=np.int32)
    y_ = np.asfortranarray(y_)
    clf = GradientBoostingClassifier(n_estimators=100, random_state=1)
    clf.fit(X, y_)
    assert_array_equal(clf.predict(T), true_result)
    assert_equal(100, len(clf.estimators_))


def test_oob_improvement():
    # Test if oob improvement has correct shape and regression test.
    clf = GradientBoostingClassifier(n_estimators=100, random_state=1,
                                     subsample=0.5)
    clf.fit(X, y)
    assert_equal(clf.oob_improvement_.shape[0], 100)
    # hard-coded regression test - change if modification in OOB computation
    assert_array_almost_equal(clf.oob_improvement_[:5],
                              np.array([0.19, 0.15, 0.12, -0.12, -0.11]),
                              decimal=2)


def test_oob_improvement_raise():
    # Test if oob improvement has correct shape.
    clf = GradientBoostingClassifier(n_estimators=100, random_state=1,
                                     subsample=1.0)
    clf.fit(X, y)
    assert_raises(AttributeError, lambda: clf.oob_improvement_)


def test_oob_multilcass_iris():
    # Check OOB improvement on multi-class dataset.
    clf = GradientBoostingClassifier(n_estimators=100, loss='deviance',
                                     random_state=1, subsample=0.5)
    clf.fit(iris.data, iris.target)
    score = clf.score(iris.data, iris.target)
    assert_greater(score, 0.9)
    assert_equal(clf.oob_improvement_.shape[0], clf.n_estimators)
    # hard-coded regression test - change if modification in OOB computation
    # FIXME: the following snippet does not yield the same results on 32 bits
    # assert_array_almost_equal(clf.oob_improvement_[:5],
    #                           np.array([12.68, 10.45, 8.18, 6.43, 5.13]),
    #                           decimal=2)


def test_verbose_output():
    # Check verbose=1 does not cause error.
    from sklearn.externals.six.moves import cStringIO as StringIO
    import sys
    old_stdout = sys.stdout
    sys.stdout = StringIO()
    clf = GradientBoostingClassifier(n_estimators=100, random_state=1,
                                     verbose=1, subsample=0.8)
    clf.fit(X, y)
    verbose_output = sys.stdout
    sys.stdout = old_stdout

    # check output
    verbose_output.seek(0)
    header = verbose_output.readline().rstrip()
    # with OOB
    true_header = ' '.join(['%10s'] + ['%16s'] * 3) % (
        'Iter', 'Train Loss', 'OOB Improve', 'Remaining Time')
    assert_equal(true_header, header)

    n_lines = sum(1 for l in verbose_output.readlines())
    # one for 1-10 and then 9 for 20-100
    assert_equal(10 + 9, n_lines)


def test_more_verbose_output():
    # Check verbose=2 does not cause error.
    from sklearn.externals.six.moves import cStringIO as StringIO
    import sys
    old_stdout = sys.stdout
    sys.stdout = StringIO()
    clf = GradientBoostingClassifier(n_estimators=100, random_state=1,
                                     verbose=2)
    clf.fit(X, y)
    verbose_output = sys.stdout
    sys.stdout = old_stdout

    # check output
    verbose_output.seek(0)
    header = verbose_output.readline().rstrip()
    # no OOB
    true_header = ' '.join(['%10s'] + ['%16s'] * 2) % (
        'Iter', 'Train Loss', 'Remaining Time')
    assert_equal(true_header, header)

    n_lines = sum(1 for l in verbose_output.readlines())
    # 100 lines for n_estimators==100
    assert_equal(100, n_lines)


def test_warm_start():
    # Test if warm start equals fit.
    X, y = datasets.make_hastie_10_2(n_samples=100, random_state=1)
    for Cls in [GradientBoostingRegressor, GradientBoostingClassifier]:
        est = Cls(n_estimators=200, max_depth=1)
        est.fit(X, y)

        est_ws = Cls(n_estimators=100, max_depth=1, warm_start=True)
        est_ws.fit(X, y)
        est_ws.set_params(n_estimators=200)
        est_ws.fit(X, y)

        assert_array_almost_equal(est_ws.predict(X), est.predict(X))


def test_warm_start_n_estimators():
    # Test if warm start equals fit - set n_estimators.
    X, y = datasets.make_hastie_10_2(n_samples=100, random_state=1)
    for Cls in [GradientBoostingRegressor, GradientBoostingClassifier]:
        est = Cls(n_estimators=300, max_depth=1)
        est.fit(X, y)

        est_ws = Cls(n_estimators=100, max_depth=1, warm_start=True)
        est_ws.fit(X, y)
        est_ws.set_params(n_estimators=300)
        est_ws.fit(X, y)

    assert_array_almost_equal(est_ws.predict(X), est.predict(X))


def test_warm_start_max_depth():
    # Test if possible to fit trees of different depth in ensemble.
    X, y = datasets.make_hastie_10_2(n_samples=100, random_state=1)
    for Cls in [GradientBoostingRegressor, GradientBoostingClassifier]:
        est = Cls(n_estimators=100, max_depth=1, warm_start=True)
        est.fit(X, y)
        est.set_params(n_estimators=110, max_depth=2)
        est.fit(X, y)

        # last 10 trees have different depth
        assert_equal(est.estimators_[0, 0].max_depth, 1)
        for i in range(1, 11):
            assert_equal(est.estimators_[-i, 0].max_depth, 2)


def test_warm_start_clear():
    # Test if fit clears state.
    X, y = datasets.make_hastie_10_2(n_samples=100, random_state=1)
    for Cls in [GradientBoostingRegressor, GradientBoostingClassifier]:
        est = Cls(n_estimators=100, max_depth=1)
        est.fit(X, y)

        est_2 = Cls(n_estimators=100, max_depth=1, warm_start=True)
        est_2.fit(X, y)  # inits state
        est_2.set_params(warm_start=False)
        est_2.fit(X, y)  # clears old state and equals est

        assert_array_almost_equal(est_2.predict(X), est.predict(X))


def test_warm_start_zero_n_estimators():
    # Test if warm start with zero n_estimators raises error
    X, y = datasets.make_hastie_10_2(n_samples=100, random_state=1)
    for Cls in [GradientBoostingRegressor, GradientBoostingClassifier]:
        est = Cls(n_estimators=100, max_depth=1, warm_start=True)
        est.fit(X, y)
        est.set_params(n_estimators=0)
        assert_raises(ValueError, est.fit, X, y)


def test_warm_start_smaller_n_estimators():
    # Test if warm start with smaller n_estimators raises error
    X, y = datasets.make_hastie_10_2(n_samples=100, random_state=1)
    for Cls in [GradientBoostingRegressor, GradientBoostingClassifier]:
        est = Cls(n_estimators=100, max_depth=1, warm_start=True)
        est.fit(X, y)
        est.set_params(n_estimators=99)
        assert_raises(ValueError, est.fit, X, y)


def test_warm_start_equal_n_estimators():
    # Test if warm start with equal n_estimators does nothing
    X, y = datasets.make_hastie_10_2(n_samples=100, random_state=1)
    for Cls in [GradientBoostingRegressor, GradientBoostingClassifier]:
        est = Cls(n_estimators=100, max_depth=1)
        est.fit(X, y)

        est2 = clone(est)
        est2.set_params(n_estimators=est.n_estimators, warm_start=True)
        est2.fit(X, y)

        assert_array_almost_equal(est2.predict(X), est.predict(X))


def test_warm_start_oob_switch():
    # Test if oob can be turned on during warm start.
    X, y = datasets.make_hastie_10_2(n_samples=100, random_state=1)
    for Cls in [GradientBoostingRegressor, GradientBoostingClassifier]:
        est = Cls(n_estimators=100, max_depth=1, warm_start=True)
        est.fit(X, y)
        est.set_params(n_estimators=110, subsample=0.5)
        est.fit(X, y)

        assert_array_equal(est.oob_improvement_[:100], np.zeros(100))
        # the last 10 are not zeros
        assert_array_equal(est.oob_improvement_[-10:] == 0.0,
                           np.zeros(10, dtype=np.bool))


def test_warm_start_oob():
    # Test if warm start OOB equals fit.
    X, y = datasets.make_hastie_10_2(n_samples=100, random_state=1)
    for Cls in [GradientBoostingRegressor, GradientBoostingClassifier]:
        est = Cls(n_estimators=200, max_depth=1, subsample=0.5,
                  random_state=1)
        est.fit(X, y)

        est_ws = Cls(n_estimators=100, max_depth=1, subsample=0.5,
                     random_state=1, warm_start=True)
        est_ws.fit(X, y)
        est_ws.set_params(n_estimators=200)
        est_ws.fit(X, y)

        assert_array_almost_equal(est_ws.oob_improvement_[:100],
                                  est.oob_improvement_[:100])


def early_stopping_monitor(i, est, locals):
    """Returns True on the 10th iteration. """
    if i == 9:
        return True
    else:
        return False


def test_monitor_early_stopping():
    # Test if monitor return value works.
    X, y = datasets.make_hastie_10_2(n_samples=100, random_state=1)

    for Cls in [GradientBoostingRegressor, GradientBoostingClassifier]:
        est = Cls(n_estimators=20, max_depth=1, random_state=1, subsample=0.5)
        est.fit(X, y, monitor=early_stopping_monitor)
        assert_equal(est.n_estimators, 20)  # this is not altered
        assert_equal(est.estimators_.shape[0], 10)
        assert_equal(est.train_score_.shape[0], 10)
        assert_equal(est.oob_improvement_.shape[0], 10)

        # try refit
        est.set_params(n_estimators=30)
        est.fit(X, y)
        assert_equal(est.n_estimators, 30)
        assert_equal(est.estimators_.shape[0], 30)
        assert_equal(est.train_score_.shape[0], 30)

        est = Cls(n_estimators=20, max_depth=1, random_state=1, subsample=0.5,
                  warm_start=True)
        est.fit(X, y, monitor=early_stopping_monitor)
        assert_equal(est.n_estimators, 20)
        assert_equal(est.estimators_.shape[0], 10)
        assert_equal(est.train_score_.shape[0], 10)
        assert_equal(est.oob_improvement_.shape[0], 10)

        # try refit
        est.set_params(n_estimators=30, warm_start=False)
        est.fit(X, y)
        assert_equal(est.n_estimators, 30)
        assert_equal(est.train_score_.shape[0], 30)
        assert_equal(est.estimators_.shape[0], 30)
        assert_equal(est.oob_improvement_.shape[0], 30)


def test_complete_classification():
    # Test greedy trees with max_depth + 1 leafs.
    from sklearn.tree._tree import TREE_LEAF
    X, y = datasets.make_hastie_10_2(n_samples=100, random_state=1)
    k = 4

    est = GradientBoostingClassifier(n_estimators=20, max_depth=None,
                                     random_state=1, max_leaf_nodes=k + 1)
    est.fit(X, y)

    tree = est.estimators_[0, 0].tree_
    assert_equal(tree.max_depth, k)
    assert_equal(tree.children_left[tree.children_left == TREE_LEAF].shape[0],
                 k + 1)


def test_complete_regression():
    # Test greedy trees with max_depth + 1 leafs.
    from sklearn.tree._tree import TREE_LEAF
    k = 4

    est = GradientBoostingRegressor(n_estimators=20, max_depth=None,
                                    random_state=1, max_leaf_nodes=k + 1)
    est.fit(boston.data, boston.target)

    tree = est.estimators_[-1, 0].tree_
    assert_equal(tree.children_left[tree.children_left == TREE_LEAF].shape[0],
                 k + 1)


def test_zero_estimator_reg():
    # Test if ZeroEstimator works for regression.
    est = GradientBoostingRegressor(n_estimators=20, max_depth=1,
                                    random_state=1, init=ZeroEstimator())
    est.fit(boston.data, boston.target)
    y_pred = est.predict(boston.data)
    mse = mean_squared_error(boston.target, y_pred)
    assert_almost_equal(mse, 33.0, decimal=0)

    est = GradientBoostingRegressor(n_estimators=20, max_depth=1,
                                    random_state=1, init='zero')
    est.fit(boston.data, boston.target)
    y_pred = est.predict(boston.data)
    mse = mean_squared_error(boston.target, y_pred)
    assert_almost_equal(mse, 33.0, decimal=0)

    est = GradientBoostingRegressor(n_estimators=20, max_depth=1,
                                    random_state=1, init='foobar')
    assert_raises(ValueError, est.fit, boston.data, boston.target)


def test_zero_estimator_clf():
    # Test if ZeroEstimator works for classification.
    X = iris.data
    y = np.array(iris.target)
    est = GradientBoostingClassifier(n_estimators=20, max_depth=1,
                                     random_state=1, init=ZeroEstimator())
    est.fit(X, y)

    assert_greater(est.score(X, y), 0.96)

    est = GradientBoostingClassifier(n_estimators=20, max_depth=1,
                                     random_state=1, init='zero')
    est.fit(X, y)

    assert_greater(est.score(X, y), 0.96)

    # binary clf
    mask = y != 0
    y[mask] = 1
    y[~mask] = 0
    est = GradientBoostingClassifier(n_estimators=20, max_depth=1,
                                     random_state=1, init='zero')
    est.fit(X, y)
    assert_greater(est.score(X, y), 0.96)

    est = GradientBoostingClassifier(n_estimators=20, max_depth=1,
                                     random_state=1, init='foobar')
    assert_raises(ValueError, est.fit, X, y)


def test_max_leaf_nodes_max_depth():
    # Test precedence of max_leaf_nodes over max_depth.
    X, y = datasets.make_hastie_10_2(n_samples=100, random_state=1)
    all_estimators = [GradientBoostingRegressor,
                      GradientBoostingClassifier]

    k = 4
    for GBEstimator in all_estimators:
        est = GBEstimator(max_depth=1, max_leaf_nodes=k).fit(X, y)
        tree = est.estimators_[0, 0].tree_
        assert_greater(tree.max_depth, 1)

        est = GBEstimator(max_depth=1).fit(X, y)
        tree = est.estimators_[0, 0].tree_
        assert_equal(tree.max_depth, 1)


def test_min_impurity_split():
    # Test if min_impurity_split of base estimators is set
    # Regression test for #8006
    X, y = datasets.make_hastie_10_2(n_samples=100, random_state=1)
    all_estimators = [GradientBoostingRegressor, GradientBoostingClassifier]

    for GBEstimator in all_estimators:
        est = GBEstimator(min_impurity_split=0.1)
        est = assert_warns_message(DeprecationWarning, "min_impurity_decrease",
                                   est.fit, X, y)
        for tree in est.estimators_.flat:
            assert_equal(tree.min_impurity_split, 0.1)


def test_min_impurity_decrease():
    X, y = datasets.make_hastie_10_2(n_samples=100, random_state=1)
    all_estimators = [GradientBoostingRegressor, GradientBoostingClassifier]

    for GBEstimator in all_estimators:
        est = GBEstimator(min_impurity_decrease=0.1)
        est.fit(X, y)
        for tree in est.estimators_.flat:
            # Simply check if the parameter is passed on correctly. Tree tests
            # will suffice for the actual working of this param
            assert_equal(tree.min_impurity_decrease, 0.1)


def test_warm_start_wo_nestimators_change():
    # Test if warm_start does nothing if n_estimators is not changed.
    # Regression test for #3513.
    clf = GradientBoostingClassifier(n_estimators=10, warm_start=True)
    clf.fit([[0, 1], [2, 3]], [0, 1])
    assert_equal(clf.estimators_.shape[0], 10)
    clf.fit([[0, 1], [2, 3]], [0, 1])
    assert_equal(clf.estimators_.shape[0], 10)


def test_probability_exponential():
    # Predict probabilities.
    clf = GradientBoostingClassifier(loss='exponential',
                                     n_estimators=100, random_state=1)

    assert_raises(ValueError, clf.predict_proba, T)

    clf.fit(X, y)
    assert_array_equal(clf.predict(T), true_result)

    # check if probabilities are in [0, 1].
    y_proba = clf.predict_proba(T)
    assert_true(np.all(y_proba >= 0.0))
    assert_true(np.all(y_proba <= 1.0))
    score = clf.decision_function(T).ravel()
    assert_array_almost_equal(y_proba[:, 1],
                              1.0 / (1.0 + np.exp(-2 * score)))

    # derive predictions from probabilities
    y_pred = clf.classes_.take(y_proba.argmax(axis=1), axis=0)
    assert_array_equal(y_pred, true_result)


def test_non_uniform_weights_toy_edge_case_reg():
    X = [[1, 0],
         [1, 0],
         [1, 0],
         [0, 1]]
    y = [0, 0, 1, 0]
    # ignore the first 2 training samples by setting their weight to 0
    sample_weight = [0, 0, 1, 1]
    for loss in ('huber', 'ls', 'lad', 'quantile'):
        gb = GradientBoostingRegressor(learning_rate=1.0, n_estimators=2,
                                       loss=loss)
        gb.fit(X, y, sample_weight=sample_weight)
        assert_greater(gb.predict([[1, 0]])[0], 0.5)


def test_non_uniform_weights_toy_edge_case_clf():
    X = [[1, 0],
         [1, 0],
         [1, 0],
         [0, 1]]
    y = [0, 0, 1, 0]
    # ignore the first 2 training samples by setting their weight to 0
    sample_weight = [0, 0, 1, 1]
    for loss in ('deviance', 'exponential'):
        gb = GradientBoostingClassifier(n_estimators=5, loss=loss)
        gb.fit(X, y, sample_weight=sample_weight)
        assert_array_equal(gb.predict([[1, 0]]), [1])


def check_sparse_input(EstimatorClass, X, X_sparse, y):
    dense = EstimatorClass(n_estimators=10, random_state=0,
                           max_depth=2).fit(X, y)
    sparse = EstimatorClass(n_estimators=10, random_state=0, max_depth=2,
                            presort=False).fit(X_sparse, y)
    auto = EstimatorClass(n_estimators=10, random_state=0, max_depth=2,
                          presort='auto').fit(X_sparse, y)

    assert_array_almost_equal(sparse.apply(X), dense.apply(X))
    assert_array_almost_equal(sparse.predict(X), dense.predict(X))
    assert_array_almost_equal(sparse.feature_importances_,
                              dense.feature_importances_)

    assert_array_almost_equal(sparse.apply(X), auto.apply(X))
    assert_array_almost_equal(sparse.predict(X), auto.predict(X))
    assert_array_almost_equal(sparse.feature_importances_,
                              auto.feature_importances_)

    assert_array_almost_equal(sparse.predict(X_sparse), dense.predict(X))
    assert_array_almost_equal(dense.predict(X_sparse), sparse.predict(X))

    if isinstance(EstimatorClass, GradientBoostingClassifier):
        assert_array_almost_equal(sparse.predict_proba(X),
                                  dense.predict_proba(X))
        assert_array_almost_equal(sparse.predict_log_proba(X),
                                  dense.predict_log_proba(X))

        assert_array_almost_equal(sparse.predict_proba(X),
                                  auto.predict_proba(X))
        assert_array_almost_equal(sparse.predict_log_proba(X),
                                  auto.predict_log_proba(X))

        assert_array_almost_equal(sparse.decision_function(X_sparse),
                                  sparse.decision_function(X))
        assert_array_almost_equal(dense.decision_function(X_sparse),
                                  sparse.decision_function(X))

        assert_array_almost_equal(
            np.array(sparse.staged_decision_function(X_sparse)),
            np.array(sparse.staged_decision_function(X)))


@skip_if_32bit
def test_sparse_input():
    ests = (GradientBoostingClassifier, GradientBoostingRegressor)
    sparse_matrices = (csr_matrix, csc_matrix, coo_matrix)

    y, X = datasets.make_multilabel_classification(random_state=0,
                                                   n_samples=50,
                                                   n_features=1,
                                                   n_classes=20)
    y = y[:, 0]

    for EstimatorClass, sparse_matrix in product(ests, sparse_matrices):
        yield check_sparse_input, EstimatorClass, X, sparse_matrix(X), y
