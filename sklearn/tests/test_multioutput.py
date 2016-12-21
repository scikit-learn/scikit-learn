import numpy as np
import scipy.sparse as sp
from sklearn.utils import shuffle
from sklearn.utils.testing import assert_almost_equal
from sklearn.utils.testing import assert_raises
from sklearn.utils.testing import assert_raises_regex
from sklearn.utils.testing import assert_array_equal
from sklearn.utils.testing import assert_equal
from sklearn.exceptions import NotFittedError
from sklearn import datasets
from sklearn.base import clone
from sklearn.ensemble import GradientBoostingRegressor, RandomForestClassifier
from sklearn.linear_model import Lasso
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
        references[:,n] = rgr.predict(X_test)

    rgr = MultiOutputRegressor(GradientBoostingRegressor(random_state=0))
    rgr.fit(X_train, y_train)
    y_pred = rgr.predict(X_test)

    assert_almost_equal(references, y_pred)


def test_multi_target_regression_one_target():
    # Test multi target regression raises
    X, y = datasets.make_regression(n_targets=1)
    X_train, y_train = X[:50], y[:50]
    X_test, y_test = X[50:], y[50:]

    rgr = MultiOutputRegressor(GradientBoostingRegressor(random_state=0))
    assert_raises(ValueError, rgr.fit, X_train, y_train)


def test_multi_target_sparse_regression():
    X, y = datasets.make_regression(n_targets=3)
    X_train, y_train = X[:50], y[:50]
    X_test, y_test = X[50:], y[50:]

    for sparse in [sp.csr_matrix, sp.csc_matrix, sp.coo_matrix, sp.dok_matrix,
                   sp.lil_matrix]:
        rgr = MultiOutputRegressor(Lasso(random_state=0))
        rgr_sparse = MultiOutputRegressor(Lasso(random_state=0))

        rgr.fit(X_train, y_train)
        rgr_sparse.fit(sparse(X_train), y_train)

        assert_almost_equal(rgr.predict(X_test), rgr_sparse.predict(sparse(X_test)))


def test_multi_target_sample_weights_api():
    X = [[1,2,3], [4,5,6]]
    y = [[3.141, 2.718], [2.718, 3.141]]
    w = [0.8, 0.6]

    rgr = MultiOutputRegressor(Lasso())
    assert_raises_regex(ValueError, "does not support sample weights",
                        rgr.fit, X, y, w)

    # no exception should be raised if the base estimator supports weights
    rgr = MultiOutputRegressor(GradientBoostingRegressor(random_state=0))
    rgr.fit(X, y, w)


def test_multi_target_sample_weights():
    # weighted regressor
    Xw = [[1,2,3], [4,5,6]]
    yw = [[3.141, 2.718], [2.718, 3.141]]
    w = [2., 1.]
    rgr_w = MultiOutputRegressor(GradientBoostingRegressor(random_state=0))
    rgr_w.fit(Xw, yw, w)

    # unweighted, but with repeated samples
    X = [[1,2,3], [1,2,3], [4,5,6]]
    y = [[3.141, 2.718], [3.141, 2.718], [2.718, 3.141]]
    rgr = MultiOutputRegressor(GradientBoostingRegressor(random_state=0))
    rgr.fit(X, y)

    X_test = [[1.5,2.5,3.5], [3.5,4.5,5.5]]
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
    assert_equal((n_samples, n_classes, n_outputs), predict_proba.shape)

    assert_array_equal(np.argmax(predict_proba, axis=1), predictions)

    # train the forest with each column and assert that predictions are equal
    for i in range(3):
        forest_ = clone(forest)  # create a clone with the same state
        forest_.fit(X, y[:, i])
        assert_equal(list(forest_.predict(X)), list(predictions[:, i]))
        assert_array_equal(list(forest_.predict_proba(X)),
                           list(predict_proba[:, :, i]))


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
