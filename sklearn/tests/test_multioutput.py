from __future__ import division
import numpy as np
import scipy.sparse as sp
from sklearn.utils import shuffle
from sklearn.utils.testing import assert_almost_equal
from sklearn.utils.testing import assert_raises
from sklearn.utils.testing import assert_false
from sklearn.utils.testing import assert_raises_regex
from sklearn.utils.testing import assert_array_equal
from sklearn.utils.testing import assert_equal
from sklearn.utils.testing import assert_not_equal
from sklearn.utils.testing import assert_array_almost_equal
from sklearn.exceptions import NotFittedError
from sklearn import datasets
from sklearn.base import clone
from sklearn.ensemble import GradientBoostingRegressor, RandomForestClassifier
from sklearn.linear_model import Lasso
from sklearn.linear_model import SGDClassifier
from sklearn.linear_model import SGDRegressor
from sklearn.linear_model import LogisticRegression
from sklearn.svm import LinearSVC
from sklearn.multiclass import OneVsRestClassifier
from sklearn.multioutput import MultiOutputRegressor, MultiOutputClassifier


def test_multi_target_regression():
    X, y = datasets.make_regression(n_targets=3)
    X_train, y_train = X[:50], y[:50]
    X_test, y_test = X[50:], y[50:]

    references = np.zeros_like(y_test)
    for n in range(3):
        rgr = GradientBoostingRegressor(random_state=0)
        rgr.fit(X_train, y_train[:, n])
        references[:, n] = rgr.predict(X_test)

    rgr = MultiOutputRegressor(GradientBoostingRegressor(random_state=0))
    rgr.fit(X_train, y_train)
    y_pred = rgr.predict(X_test)

    assert_almost_equal(references, y_pred)


def test_multi_target_regression_partial_fit():
    X, y = datasets.make_regression(n_targets=3)
    X_train, y_train = X[:50], y[:50]
    X_test, y_test = X[50:], y[50:]

    references = np.zeros_like(y_test)
    half_index = 25
    for n in range(3):
        sgr = SGDRegressor(random_state=0)
        sgr.partial_fit(X_train[:half_index], y_train[:half_index, n])
        sgr.partial_fit(X_train[half_index:], y_train[half_index:, n])
        references[:, n] = sgr.predict(X_test)

    sgr = MultiOutputRegressor(SGDRegressor(random_state=0))

    sgr.partial_fit(X_train[:half_index], y_train[:half_index])
    sgr.partial_fit(X_train[half_index:], y_train[half_index:])

    y_pred = sgr.predict(X_test)
    assert_almost_equal(references, y_pred)


def test_multi_target_regression_one_target():
    # Test multi target regression raises
    X, y = datasets.make_regression(n_targets=1)

    rgr = MultiOutputRegressor(GradientBoostingRegressor(random_state=0))
    assert_raises(ValueError, rgr.fit, X, y)


def test_multi_target_sparse_regression():
    X, y = datasets.make_regression(n_targets=3)
    X_train, y_train = X[:50], y[:50]
    X_test = X[50:]

    for sparse in [sp.csr_matrix, sp.csc_matrix, sp.coo_matrix, sp.dok_matrix,
                   sp.lil_matrix]:
        rgr = MultiOutputRegressor(Lasso(random_state=0))
        rgr_sparse = MultiOutputRegressor(Lasso(random_state=0))

        rgr.fit(X_train, y_train)
        rgr_sparse.fit(sparse(X_train), y_train)

        assert_almost_equal(rgr.predict(X_test),
                            rgr_sparse.predict(sparse(X_test)))


def test_multi_target_sample_weights_api():
    X = [[1, 2, 3], [4, 5, 6]]
    y = [[3.141, 2.718], [2.718, 3.141]]
    w = [0.8, 0.6]

    rgr = MultiOutputRegressor(Lasso())
    assert_raises_regex(ValueError, "does not support sample weights",
                        rgr.fit, X, y, w)

    # no exception should be raised if the base estimator supports weights
    rgr = MultiOutputRegressor(GradientBoostingRegressor(random_state=0))
    rgr.fit(X, y, w)


def test_multi_target_sample_weight_partial_fit():
    # weighted regressor
    X = [[1, 2, 3], [4, 5, 6]]
    y = [[3.141, 2.718], [2.718, 3.141]]
    w = [2., 1.]
    rgr_w = MultiOutputRegressor(SGDRegressor(random_state=0))
    rgr_w.partial_fit(X, y, w)

    # weighted with different weights
    w = [2., 2.]
    rgr = MultiOutputRegressor(SGDRegressor(random_state=0))
    rgr.partial_fit(X, y, w)

    assert_not_equal(rgr.predict(X)[0][0], rgr_w.predict(X)[0][0])


def test_multi_target_sample_weights():
    # weighted regressor
    Xw = [[1, 2, 3], [4, 5, 6]]
    yw = [[3.141, 2.718], [2.718, 3.141]]
    w = [2., 1.]
    rgr_w = MultiOutputRegressor(GradientBoostingRegressor(random_state=0))
    rgr_w.fit(Xw, yw, w)

    # unweighted, but with repeated samples
    X = [[1, 2, 3], [1, 2, 3], [4, 5, 6]]
    y = [[3.141, 2.718], [3.141, 2.718], [2.718, 3.141]]
    rgr = MultiOutputRegressor(GradientBoostingRegressor(random_state=0))
    rgr.fit(X, y)

    X_test = [[1.5, 2.5, 3.5], [3.5, 4.5, 5.5]]
    assert_almost_equal(rgr.predict(X_test), rgr_w.predict(X_test))


# Import the data
iris = datasets.load_iris()
# create a multiple targets by randomized shuffling and concatenating y.
X = iris.data
y1 = iris.target
y2 = shuffle(y1, random_state=1)
y3 = shuffle(y1, random_state=2)
y = np.column_stack((y1, y2, y3))
n_samples, n_features = X.shape
n_outputs = y.shape[1]
n_classes = len(np.unique(y1))
classes = list(map(np.unique, (y1, y2, y3)))


def test_multi_output_classification_partial_fit_parallelism():
    sgd_linear_clf = SGDClassifier(loss='log', random_state=1)
    mor = MultiOutputClassifier(sgd_linear_clf, n_jobs=-1)
    mor.partial_fit(X, y, classes)
    est1 = mor.estimators_[0]
    mor.partial_fit(X, y)
    est2 = mor.estimators_[0]
    # parallelism requires this to be the case for a sane implementation
    assert_false(est1 is est2)


def test_multi_output_classification_partial_fit():
    # test if multi_target initializes correctly with base estimator and fit
    # assert predictions work as expected for predict

    sgd_linear_clf = SGDClassifier(loss='log', random_state=1)
    multi_target_linear = MultiOutputClassifier(sgd_linear_clf)

    # train the multi_target_linear and also get the predictions.
    half_index = X.shape[0] // 2
    multi_target_linear.partial_fit(
        X[:half_index], y[:half_index], classes=classes)

    first_predictions = multi_target_linear.predict(X)
    assert_equal((n_samples, n_outputs), first_predictions.shape)

    multi_target_linear.partial_fit(X[half_index:], y[half_index:])
    second_predictions = multi_target_linear.predict(X)
    assert_equal((n_samples, n_outputs), second_predictions.shape)

    # train the linear classification with each column and assert that
    # predictions are equal after first partial_fit and second partial_fit
    for i in range(3):
        # create a clone with the same state
        sgd_linear_clf = clone(sgd_linear_clf)
        sgd_linear_clf.partial_fit(
            X[:half_index], y[:half_index, i], classes=classes[i])
        assert_array_equal(sgd_linear_clf.predict(X), first_predictions[:, i])
        sgd_linear_clf.partial_fit(X[half_index:], y[half_index:, i])
        assert_array_equal(sgd_linear_clf.predict(X), second_predictions[:, i])


def test_mutli_output_classifiation_partial_fit_no_first_classes_exception():
    sgd_linear_clf = SGDClassifier(loss='log', random_state=1)
    multi_target_linear = MultiOutputClassifier(sgd_linear_clf)
    assert_raises_regex(ValueError, "classes must be passed on the first call "
                                    "to partial_fit.",
                        multi_target_linear.partial_fit, X, y)


def test_multi_output_classification():
    # test if multi_target initializes correctly with base estimator and fit
    # assert predictions work as expected for predict, prodict_proba and score

    forest = RandomForestClassifier(n_estimators=10, random_state=1)
    multi_target_forest = MultiOutputClassifier(forest)

    # train the multi_target_forest and also get the predictions.
    multi_target_forest.fit(X, y)

    predictions = multi_target_forest.predict(X)
    assert_equal((n_samples, n_outputs), predictions.shape)

    predict_proba = multi_target_forest.predict_proba(X)

    assert len(predict_proba) == n_outputs
    for class_probabilities in predict_proba:
        assert_equal((n_samples, n_classes), class_probabilities.shape)

    assert_array_equal(np.argmax(np.dstack(predict_proba), axis=1),
                       predictions)

    # train the forest with each column and assert that predictions are equal
    for i in range(3):
        forest_ = clone(forest)  # create a clone with the same state
        forest_.fit(X, y[:, i])
        assert_equal(list(forest_.predict(X)), list(predictions[:, i]))
        assert_array_equal(list(forest_.predict_proba(X)),
                           list(predict_proba[i]))


def test_multiclass_multioutput_estimator():
    # test to check meta of meta estimators
    svc = LinearSVC(random_state=0)
    multi_class_svc = OneVsRestClassifier(svc)
    multi_target_svc = MultiOutputClassifier(multi_class_svc)

    multi_target_svc.fit(X, y)

    predictions = multi_target_svc.predict(X)
    assert_equal((n_samples, n_outputs), predictions.shape)

    # train the forest with each column and assert that predictions are equal
    for i in range(3):
        multi_class_svc_ = clone(multi_class_svc)  # create a clone
        multi_class_svc_.fit(X, y[:, i])
        assert_equal(list(multi_class_svc_.predict(X)),
                     list(predictions[:, i]))


def test_multiclass_multioutput_estimator_predict_proba():
    seed = 542

    # make test deterministic
    rng = np.random.RandomState(seed)

    # random features
    X = rng.normal(size=(5, 5))

    # random labels
    y1 = np.array(['b', 'a', 'a', 'b', 'a']).reshape(5, 1)  # 2 classes
    y2 = np.array(['d', 'e', 'f', 'e', 'd']).reshape(5, 1)  # 3 classes

    Y = np.concatenate([y1, y2], axis=1)

    clf = MultiOutputClassifier(LogisticRegression(random_state=seed))

    clf.fit(X, Y)

    y_result = clf.predict_proba(X)
    y_actual = [np.array([[0.23481764, 0.76518236],
                          [0.67196072, 0.32803928],
                          [0.54681448, 0.45318552],
                          [0.34883923, 0.65116077],
                          [0.73687069, 0.26312931]]),
                np.array([[0.5171785, 0.23878628, 0.24403522],
                          [0.22141451, 0.64102704, 0.13755846],
                          [0.16751315, 0.18256843, 0.64991843],
                          [0.27357372, 0.55201592, 0.17441036],
                          [0.65745193, 0.26062899, 0.08191907]])]

    for i in range(len(y_actual)):
        assert_almost_equal(y_result[i], y_actual[i])


def test_multi_output_classification_sample_weights():
    # weighted classifier
    Xw = [[1, 2, 3], [4, 5, 6]]
    yw = [[3, 2], [2, 3]]
    w = np.asarray([2., 1.])
    forest = RandomForestClassifier(n_estimators=10, random_state=1)
    clf_w = MultiOutputClassifier(forest)
    clf_w.fit(Xw, yw, w)

    # unweighted, but with repeated samples
    X = [[1, 2, 3], [1, 2, 3], [4, 5, 6]]
    y = [[3, 2], [3, 2], [2, 3]]
    forest = RandomForestClassifier(n_estimators=10, random_state=1)
    clf = MultiOutputClassifier(forest)
    clf.fit(X, y)

    X_test = [[1.5, 2.5, 3.5], [3.5, 4.5, 5.5]]
    assert_almost_equal(clf.predict(X_test), clf_w.predict(X_test))


def test_multi_output_classification_partial_fit_sample_weights():
    # weighted classifier
    Xw = [[1, 2, 3], [4, 5, 6], [1.5, 2.5, 3.5]]
    yw = [[3, 2], [2, 3], [3, 2]]
    w = np.asarray([2., 1., 1.])
    sgd_linear_clf = SGDClassifier(random_state=1)
    clf_w = MultiOutputClassifier(sgd_linear_clf)
    clf_w.fit(Xw, yw, w)

    # unweighted, but with repeated samples
    X = [[1, 2, 3], [1, 2, 3], [4, 5, 6], [1.5, 2.5, 3.5]]
    y = [[3, 2], [3, 2], [2, 3], [3, 2]]
    sgd_linear_clf = SGDClassifier(random_state=1)
    clf = MultiOutputClassifier(sgd_linear_clf)
    clf.fit(X, y)
    X_test = [[1.5, 2.5, 3.5]]
    assert_array_almost_equal(clf.predict(X_test), clf_w.predict(X_test))


def test_multi_output_exceptions():
    # NotFittedError when fit is not done but score, predict and
    # and predict_proba are called
    moc = MultiOutputClassifier(LinearSVC(random_state=0))
    assert_raises(NotFittedError, moc.predict, y)
    assert_raises(NotFittedError, moc.predict_proba, y)
    assert_raises(NotFittedError, moc.score, X, y)
    # ValueError when number of outputs is different
    # for fit and score
    y_new = np.column_stack((y1, y2))
    moc.fit(X, y)
    assert_raises(ValueError, moc.score, X, y_new)
