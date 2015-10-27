"""
Testing for Support Vector Machine module (sklearn.svm)

TODO: remove hard coded numerical results when possible
"""

import numpy as np
import itertools

from numpy.testing import assert_array_equal, assert_array_almost_equal
from numpy.testing import assert_almost_equal
from scipy import sparse
from nose.tools import assert_raises, assert_true, assert_equal, assert_false

from sklearn.base import ChangedBehaviorWarning
from sklearn import svm, linear_model, datasets, metrics, base
from sklearn.cross_validation import train_test_split
from sklearn.datasets import make_classification, make_blobs
from sklearn.metrics import f1_score
from sklearn.metrics.pairwise import rbf_kernel
from sklearn.utils import check_random_state
from sklearn.utils import ConvergenceWarning
from sklearn.utils.validation import NotFittedError
from sklearn.utils.testing import assert_greater, assert_in, assert_less
from sklearn.utils.testing import assert_raises_regexp, assert_warns
from sklearn.utils.testing import assert_warns_message, assert_raise_message
from sklearn.utils.testing import ignore_warnings

from sklearn.multiclass import OneVsRestClassifier

# toy sample
X = [[-2, -1], [-1, -1], [-1, -2], [1, 1], [1, 2], [2, 1]]
Y = [1, 1, 1, 2, 2, 2]
T = [[-1, -1], [2, 2], [3, 2]]
true_result = [1, 2, 2]

# also load the iris dataset
iris = datasets.load_iris()
rng = check_random_state(42)
perm = rng.permutation(iris.target.size)
iris.data = iris.data[perm]
iris.target = iris.target[perm]


def test_libsvm_parameters():
    # Test parameters on classes that make use of libsvm.
    clf = svm.SVC(kernel='linear').fit(X, Y)
    assert_array_equal(clf.dual_coef_, [[-0.25, .25]])
    assert_array_equal(clf.support_, [1, 3])
    assert_array_equal(clf.support_vectors_, (X[1], X[3]))
    assert_array_equal(clf.intercept_, [0.])
    assert_array_equal(clf.predict(X), Y)


def test_libsvm_iris():
    # Check consistency on dataset iris.

    # shuffle the dataset so that labels are not ordered
    for k in ('linear', 'rbf'):
        clf = svm.SVC(kernel=k).fit(iris.data, iris.target)
        assert_greater(np.mean(clf.predict(iris.data) == iris.target), 0.9)

    assert_array_equal(clf.classes_, np.sort(clf.classes_))

    # check also the low-level API
    model = svm.libsvm.fit(iris.data, iris.target.astype(np.float64))
    pred = svm.libsvm.predict(iris.data, *model)
    assert_greater(np.mean(pred == iris.target), .95)

    model = svm.libsvm.fit(iris.data, iris.target.astype(np.float64),
                           kernel='linear')
    pred = svm.libsvm.predict(iris.data, *model, kernel='linear')
    assert_greater(np.mean(pred == iris.target), .95)

    pred = svm.libsvm.cross_validation(iris.data,
                                       iris.target.astype(np.float64), 5,
                                       kernel='linear',
                                       random_seed=0)
    assert_greater(np.mean(pred == iris.target), .95)

    # If random_seed >= 0, the libsvm rng is seeded (by calling `srand`), hence
    # we should get deteriministic results (assuming that there is no other
    # thread calling this wrapper calling `srand` concurrently).
    pred2 = svm.libsvm.cross_validation(iris.data,
                                        iris.target.astype(np.float64), 5,
                                        kernel='linear',
                                        random_seed=0)
    assert_array_equal(pred, pred2)


@ignore_warnings
def test_single_sample_1d():
    # Test whether SVCs work on a single sample given as a 1-d array

    clf = svm.SVC().fit(X, Y)
    clf.predict(X[0])

    clf = svm.LinearSVC(random_state=0).fit(X, Y)
    clf.predict(X[0])


def test_precomputed():
    # SVC with a precomputed kernel.
    # We test it with a toy dataset and with iris.
    clf = svm.SVC(kernel='precomputed')
    # Gram matrix for train data (square matrix)
    # (we use just a linear kernel)
    K = np.dot(X, np.array(X).T)
    clf.fit(K, Y)
    # Gram matrix for test data (rectangular matrix)
    KT = np.dot(T, np.array(X).T)
    pred = clf.predict(KT)
    assert_raises(ValueError, clf.predict, KT.T)

    assert_array_equal(clf.dual_coef_, [[-0.25, .25]])
    assert_array_equal(clf.support_, [1, 3])
    assert_array_equal(clf.intercept_, [0])
    assert_array_almost_equal(clf.support_, [1, 3])
    assert_array_equal(pred, true_result)

    # Gram matrix for test data but compute KT[i,j]
    # for support vectors j only.
    KT = np.zeros_like(KT)
    for i in range(len(T)):
        for j in clf.support_:
            KT[i, j] = np.dot(T[i], X[j])

    pred = clf.predict(KT)
    assert_array_equal(pred, true_result)

    # same as before, but using a callable function instead of the kernel
    # matrix. kernel is just a linear kernel

    kfunc = lambda x, y: np.dot(x, y.T)
    clf = svm.SVC(kernel=kfunc)
    clf.fit(X, Y)
    pred = clf.predict(T)

    assert_array_equal(clf.dual_coef_, [[-0.25, .25]])
    assert_array_equal(clf.intercept_, [0])
    assert_array_almost_equal(clf.support_, [1, 3])
    assert_array_equal(pred, true_result)

    # test a precomputed kernel with the iris dataset
    # and check parameters against a linear SVC
    clf = svm.SVC(kernel='precomputed')
    clf2 = svm.SVC(kernel='linear')
    K = np.dot(iris.data, iris.data.T)
    clf.fit(K, iris.target)
    clf2.fit(iris.data, iris.target)
    pred = clf.predict(K)
    assert_array_almost_equal(clf.support_, clf2.support_)
    assert_array_almost_equal(clf.dual_coef_, clf2.dual_coef_)
    assert_array_almost_equal(clf.intercept_, clf2.intercept_)
    assert_almost_equal(np.mean(pred == iris.target), .99, decimal=2)

    # Gram matrix for test data but compute KT[i,j]
    # for support vectors j only.
    K = np.zeros_like(K)
    for i in range(len(iris.data)):
        for j in clf.support_:
            K[i, j] = np.dot(iris.data[i], iris.data[j])

    pred = clf.predict(K)
    assert_almost_equal(np.mean(pred == iris.target), .99, decimal=2)

    clf = svm.SVC(kernel=kfunc)
    clf.fit(iris.data, iris.target)
    assert_almost_equal(np.mean(pred == iris.target), .99, decimal=2)


def test_svr():
    # Test Support Vector Regression

    diabetes = datasets.load_diabetes()
    for clf in (svm.NuSVR(kernel='linear', nu=.4, C=1.0),
                svm.NuSVR(kernel='linear', nu=.4, C=10.),
                svm.SVR(kernel='linear', C=10.),
                svm.LinearSVR(C=10.),
                svm.LinearSVR(C=10.),
                ):
        clf.fit(diabetes.data, diabetes.target)
        assert_greater(clf.score(diabetes.data, diabetes.target), 0.02)

    # non-regression test; previously, BaseLibSVM would check that
    # len(np.unique(y)) < 2, which must only be done for SVC
    svm.SVR().fit(diabetes.data, np.ones(len(diabetes.data)))
    svm.LinearSVR().fit(diabetes.data, np.ones(len(diabetes.data)))


def test_linearsvr():
    # check that SVR(kernel='linear') and LinearSVC() give
    # comparable results
    diabetes = datasets.load_diabetes()
    lsvr = svm.LinearSVR(C=1e3).fit(diabetes.data, diabetes.target)
    score1 = lsvr.score(diabetes.data, diabetes.target)

    svr = svm.SVR(kernel='linear', C=1e3).fit(diabetes.data, diabetes.target)
    score2 = svr.score(diabetes.data, diabetes.target)

    assert np.linalg.norm(lsvr.coef_ - svr.coef_) / np.linalg.norm(svr.coef_) < .1
    assert np.abs(score1 - score2) < 0.1


def test_svr_errors():
    X = [[0.0], [1.0]]
    y = [0.0, 0.5]

    # Bad kernel
    clf = svm.SVR(kernel=lambda x, y: np.array([[1.0]]))
    clf.fit(X, y)
    assert_raises(ValueError, clf.predict, X)


def test_oneclass():
    # Test OneClassSVM
    clf = svm.OneClassSVM()
    clf.fit(X)
    pred = clf.predict(T)

    assert_array_almost_equal(pred, [-1, -1, -1])
    assert_array_almost_equal(clf.intercept_, [-1.008], decimal=3)
    assert_array_almost_equal(clf.dual_coef_,
                              [[0.632, 0.233, 0.633, 0.234, 0.632, 0.633]],
                              decimal=3)
    assert_raises(ValueError, lambda: clf.coef_)


def test_oneclass_decision_function():
    # Test OneClassSVM decision function
    clf = svm.OneClassSVM()
    rnd = check_random_state(2)

    # Generate train data
    X = 0.3 * rnd.randn(100, 2)
    X_train = np.r_[X + 2, X - 2]

    # Generate some regular novel observations
    X = 0.3 * rnd.randn(20, 2)
    X_test = np.r_[X + 2, X - 2]
    # Generate some abnormal novel observations
    X_outliers = rnd.uniform(low=-4, high=4, size=(20, 2))

    # fit the model
    clf = svm.OneClassSVM(nu=0.1, kernel="rbf", gamma=0.1)
    clf.fit(X_train)

    # predict things
    y_pred_test = clf.predict(X_test)
    assert_greater(np.mean(y_pred_test == 1), .9)
    y_pred_outliers = clf.predict(X_outliers)
    assert_greater(np.mean(y_pred_outliers == -1), .9)
    dec_func_test = clf.decision_function(X_test)
    assert_array_equal((dec_func_test > 0).ravel(), y_pred_test == 1)
    dec_func_outliers = clf.decision_function(X_outliers)
    assert_array_equal((dec_func_outliers > 0).ravel(), y_pred_outliers == 1)


def test_tweak_params():
    # Make sure some tweaking of parameters works.
    # We change clf.dual_coef_ at run time and expect .predict() to change
    # accordingly. Notice that this is not trivial since it involves a lot
    # of C/Python copying in the libsvm bindings.
    # The success of this test ensures that the mapping between libsvm and
    # the python classifier is complete.
    clf = svm.SVC(kernel='linear', C=1.0)
    clf.fit(X, Y)
    assert_array_equal(clf.dual_coef_, [[-.25, .25]])
    assert_array_equal(clf.predict([[-.1, -.1]]), [1])
    clf._dual_coef_ = np.array([[.0, 1.]])
    assert_array_equal(clf.predict([[-.1, -.1]]), [2])


def test_probability():
    # Predict probabilities using SVC
    # This uses cross validation, so we use a slightly bigger testing set.

    for clf in (svm.SVC(probability=True, random_state=0, C=1.0),
                svm.NuSVC(probability=True, random_state=0)):

        clf.fit(iris.data, iris.target)

        prob_predict = clf.predict_proba(iris.data)
        assert_array_almost_equal(
            np.sum(prob_predict, 1), np.ones(iris.data.shape[0]))
        assert_true(np.mean(np.argmax(prob_predict, 1)
                    == clf.predict(iris.data)) > 0.9)

        assert_almost_equal(clf.predict_proba(iris.data),
                            np.exp(clf.predict_log_proba(iris.data)), 8)


def test_decision_function():
    # Test decision_function
    # Sanity check, test that decision_function implemented in python
    # returns the same as the one in libsvm
    # multi class:
    clf = svm.SVC(kernel='linear', C=0.1,
                  decision_function_shape='ovo').fit(iris.data, iris.target)

    dec = np.dot(iris.data, clf.coef_.T) + clf.intercept_

    assert_array_almost_equal(dec, clf.decision_function(iris.data))

    # binary:
    clf.fit(X, Y)
    dec = np.dot(X, clf.coef_.T) + clf.intercept_
    prediction = clf.predict(X)
    assert_array_almost_equal(dec.ravel(), clf.decision_function(X))
    assert_array_almost_equal(
        prediction,
        clf.classes_[(clf.decision_function(X) > 0).astype(np.int)])
    expected = np.array([-1., -0.66, -1., 0.66, 1., 1.])
    assert_array_almost_equal(clf.decision_function(X), expected, 2)

    # kernel binary:
    clf = svm.SVC(kernel='rbf', gamma=1, decision_function_shape='ovo')
    clf.fit(X, Y)

    rbfs = rbf_kernel(X, clf.support_vectors_, gamma=clf.gamma)
    dec = np.dot(rbfs, clf.dual_coef_.T) + clf.intercept_
    assert_array_almost_equal(dec.ravel(), clf.decision_function(X))


def test_decision_function_shape():
    # check that decision_function_shape='ovr' gives
    # correct shape and is consistent with predict

    clf = svm.SVC(kernel='linear', C=0.1,
                  decision_function_shape='ovr').fit(iris.data, iris.target)
    dec = clf.decision_function(iris.data)
    assert_equal(dec.shape, (len(iris.data), 3))
    assert_array_equal(clf.predict(iris.data), np.argmax(dec, axis=1))

    # with five classes:
    X, y = make_blobs(n_samples=80, centers=5, random_state=0)
    X_train, X_test, y_train, y_test = train_test_split(X, y, random_state=0)

    clf = svm.SVC(kernel='linear', C=0.1,
                  decision_function_shape='ovr').fit(X_train, y_train)
    dec = clf.decision_function(X_test)
    assert_equal(dec.shape, (len(X_test), 5))
    assert_array_equal(clf.predict(X_test), np.argmax(dec, axis=1))

    # check shape of ovo_decition_function=True
    clf = svm.SVC(kernel='linear', C=0.1,
                  decision_function_shape='ovo').fit(X_train, y_train)
    dec = clf.decision_function(X_train)
    assert_equal(dec.shape, (len(X_train), 10))

    # check deprecation warning
    clf = svm.SVC(kernel='linear', C=0.1).fit(X_train, y_train)
    msg = "change the shape of the decision function"
    dec = assert_warns_message(ChangedBehaviorWarning, msg,
                               clf.decision_function, X_train)
    assert_equal(dec.shape, (len(X_train), 10))


def test_svr_predict():
    # Test SVR's decision_function
    # Sanity check, test that predict implemented in python
    # returns the same as the one in libsvm

    X = iris.data
    y = iris.target

    # linear kernel
    reg = svm.SVR(kernel='linear', C=0.1).fit(X, y)

    dec = np.dot(X, reg.coef_.T) + reg.intercept_
    assert_array_almost_equal(dec.ravel(), reg.predict(X).ravel())

    # rbf kernel
    reg = svm.SVR(kernel='rbf', gamma=1).fit(X, y)

    rbfs = rbf_kernel(X, reg.support_vectors_, gamma=reg.gamma)
    dec = np.dot(rbfs, reg.dual_coef_.T) + reg.intercept_
    assert_array_almost_equal(dec.ravel(), reg.predict(X).ravel())


def test_weight():
    # Test class weights
    clf = svm.SVC(class_weight={1: 0.1})
    # we give a small weights to class 1
    clf.fit(X, Y)
    # so all predicted values belong to class 2
    assert_array_almost_equal(clf.predict(X), [2] * 6)

    X_, y_ = make_classification(n_samples=200, n_features=10,
                                 weights=[0.833, 0.167], random_state=2)

    for clf in (linear_model.LogisticRegression(),
                svm.LinearSVC(random_state=0), svm.SVC()):
        clf.set_params(class_weight={0: .1, 1: 10})
        clf.fit(X_[:100], y_[:100])
        y_pred = clf.predict(X_[100:])
        assert_true(f1_score(y_[100:], y_pred) > .3)


def test_sample_weights():
    # Test weights on individual samples
    # TODO: check on NuSVR, OneClass, etc.
    clf = svm.SVC()
    clf.fit(X, Y)
    assert_array_equal(clf.predict([X[2]]), [1.])

    sample_weight = [.1] * 3 + [10] * 3
    clf.fit(X, Y, sample_weight=sample_weight)
    assert_array_equal(clf.predict([X[2]]), [2.])

    # test that rescaling all samples is the same as changing C
    clf = svm.SVC()
    clf.fit(X, Y)
    dual_coef_no_weight = clf.dual_coef_
    clf.set_params(C=100)
    clf.fit(X, Y, sample_weight=np.repeat(0.01, len(X)))
    assert_array_almost_equal(dual_coef_no_weight, clf.dual_coef_)


def test_auto_weight():
    # Test class weights for imbalanced data
    from sklearn.linear_model import LogisticRegression
    # We take as dataset the two-dimensional projection of iris so
    # that it is not separable and remove half of predictors from
    # class 1.
    # We add one to the targets as a non-regression test: class_weight="balanced"
    # used to work only when the labels where a range [0..K).
    from sklearn.utils import compute_class_weight
    X, y = iris.data[:, :2], iris.target + 1
    unbalanced = np.delete(np.arange(y.size), np.where(y > 2)[0][::2])

    classes = np.unique(y[unbalanced])
    class_weights = compute_class_weight('balanced', classes, y[unbalanced])
    assert_true(np.argmax(class_weights) == 2)

    for clf in (svm.SVC(kernel='linear'), svm.LinearSVC(random_state=0),
                LogisticRegression()):
        # check that score is better when class='balanced' is set.
        y_pred = clf.fit(X[unbalanced], y[unbalanced]).predict(X)
        clf.set_params(class_weight='balanced')
        y_pred_balanced = clf.fit(X[unbalanced], y[unbalanced],).predict(X)
        assert_true(metrics.f1_score(y, y_pred, average='weighted')
                    <= metrics.f1_score(y, y_pred_balanced,
                                        average='weighted'))


def test_bad_input():
    # Test that it gives proper exception on deficient input
    # impossible value of C
    assert_raises(ValueError, svm.SVC(C=-1).fit, X, Y)

    # impossible value of nu
    clf = svm.NuSVC(nu=0.0)
    assert_raises(ValueError, clf.fit, X, Y)

    Y2 = Y[:-1]  # wrong dimensions for labels
    assert_raises(ValueError, clf.fit, X, Y2)

    # Test with arrays that are non-contiguous.
    for clf in (svm.SVC(), svm.LinearSVC(random_state=0)):
        Xf = np.asfortranarray(X)
        assert_false(Xf.flags['C_CONTIGUOUS'])
        yf = np.ascontiguousarray(np.tile(Y, (2, 1)).T)
        yf = yf[:, -1]
        assert_false(yf.flags['F_CONTIGUOUS'])
        assert_false(yf.flags['C_CONTIGUOUS'])
        clf.fit(Xf, yf)
        assert_array_equal(clf.predict(T), true_result)

    # error for precomputed kernelsx
    clf = svm.SVC(kernel='precomputed')
    assert_raises(ValueError, clf.fit, X, Y)

    # sample_weight bad dimensions
    clf = svm.SVC()
    assert_raises(ValueError, clf.fit, X, Y, sample_weight=range(len(X) - 1))

    # predict with sparse input when trained with dense
    clf = svm.SVC().fit(X, Y)
    assert_raises(ValueError, clf.predict, sparse.lil_matrix(X))

    Xt = np.array(X).T
    clf.fit(np.dot(X, Xt), Y)
    assert_raises(ValueError, clf.predict, X)

    clf = svm.SVC()
    clf.fit(X, Y)
    assert_raises(ValueError, clf.predict, Xt)


def test_sparse_precomputed():
    clf = svm.SVC(kernel='precomputed')
    sparse_gram = sparse.csr_matrix([[1, 0], [0, 1]])
    try:
        clf.fit(sparse_gram, [0, 1])
        assert not "reached"
    except TypeError as e:
        assert_in("Sparse precomputed", str(e))


def test_linearsvc_parameters():
    # Test possible parameter combinations in LinearSVC
    # Generate list of possible parameter combinations
    losses = ['hinge', 'squared_hinge', 'logistic_regression', 'foo']
    penalties, duals = ['l1', 'l2', 'bar'], [True, False]

    X, y = make_classification(n_samples=5, n_features=5)

    for loss, penalty, dual in itertools.product(losses, penalties, duals):
        clf = svm.LinearSVC(penalty=penalty, loss=loss, dual=dual)
        if ((loss, penalty) == ('hinge', 'l1') or
                (loss, penalty, dual) == ('hinge', 'l2', False) or
                (penalty, dual) == ('l1', True) or
                loss == 'foo' or penalty == 'bar'):

            assert_raises_regexp(ValueError,
                                 "Unsupported set of arguments.*penalty='%s.*"
                                 "loss='%s.*dual=%s"
                                 % (penalty, loss, dual),
                                 clf.fit, X, y)
        else:
            clf.fit(X, y)

    # Incorrect loss value - test if explicit error message is raised
    assert_raises_regexp(ValueError, ".*loss='l3' is not supported.*",
                         svm.LinearSVC(loss="l3").fit, X, y)


# FIXME remove in 1.0
def test_linearsvx_loss_penalty_deprecations():
    X, y = [[0.0], [1.0]], [0, 1]

    msg = ("loss='%s' has been deprecated in favor of "
           "loss='%s' as of 0.16. Backward compatibility"
           " for the %s will be removed in %s")

    # LinearSVC
    # loss l1/L1 --> hinge
    assert_warns_message(DeprecationWarning,
                         msg % ("l1", "hinge", "loss='l1'", "1.0"),
                         svm.LinearSVC(loss="l1").fit, X, y)

    # loss l2/L2 --> squared_hinge
    assert_warns_message(DeprecationWarning,
                         msg % ("L2", "squared_hinge", "loss='L2'", "1.0"),
                         svm.LinearSVC(loss="L2").fit, X, y)

    # LinearSVR
    # loss l1/L1 --> epsilon_insensitive
    assert_warns_message(DeprecationWarning,
                         msg % ("L1", "epsilon_insensitive", "loss='L1'",
                                "1.0"),
                         svm.LinearSVR(loss="L1").fit, X, y)

    # loss l2/L2 --> squared_epsilon_insensitive
    assert_warns_message(DeprecationWarning,
                         msg % ("l2", "squared_epsilon_insensitive",
                                "loss='l2'", "1.0"),
                         svm.LinearSVR(loss="l2").fit, X, y)


# FIXME remove in 0.18
def test_linear_svx_uppercase_loss_penalty():
    # Check if Upper case notation is supported by _fit_liblinear
    # which is called by fit
    X, y = [[0.0], [1.0]], [0, 1]

    msg = ("loss='%s' has been deprecated in favor of "
           "loss='%s' as of 0.16. Backward compatibility"
           " for the uppercase notation will be removed in %s")

    # loss SQUARED_hinge --> squared_hinge
    assert_warns_message(DeprecationWarning,
                         msg % ("SQUARED_hinge", "squared_hinge", "0.18"),
                         svm.LinearSVC(loss="SQUARED_hinge").fit, X, y)

    # penalty L2 --> l2
    assert_warns_message(DeprecationWarning,
                         msg.replace("loss", "penalty")
                         % ("L2", "l2", "0.18"),
                         svm.LinearSVC(penalty="L2").fit, X, y)

    # loss EPSILON_INSENSITIVE --> epsilon_insensitive
    assert_warns_message(DeprecationWarning,
                         msg % ("EPSILON_INSENSITIVE", "epsilon_insensitive",
                                "0.18"),
                         svm.LinearSVR(loss="EPSILON_INSENSITIVE").fit, X, y)


def test_linearsvc():
    # Test basic routines using LinearSVC
    clf = svm.LinearSVC(random_state=0).fit(X, Y)

    # by default should have intercept
    assert_true(clf.fit_intercept)

    assert_array_equal(clf.predict(T), true_result)
    assert_array_almost_equal(clf.intercept_, [0], decimal=3)

    # the same with l1 penalty
    clf = svm.LinearSVC(penalty='l1', loss='squared_hinge', dual=False, random_state=0).fit(X, Y)
    assert_array_equal(clf.predict(T), true_result)

    # l2 penalty with dual formulation
    clf = svm.LinearSVC(penalty='l2', dual=True, random_state=0).fit(X, Y)
    assert_array_equal(clf.predict(T), true_result)

    # l2 penalty, l1 loss
    clf = svm.LinearSVC(penalty='l2', loss='hinge', dual=True, random_state=0)
    clf.fit(X, Y)
    assert_array_equal(clf.predict(T), true_result)

    # test also decision function
    dec = clf.decision_function(T)
    res = (dec > 0).astype(np.int) + 1
    assert_array_equal(res, true_result)


def test_linearsvc_crammer_singer():
    # Test LinearSVC with crammer_singer multi-class svm
    ovr_clf = svm.LinearSVC(random_state=0).fit(iris.data, iris.target)
    cs_clf = svm.LinearSVC(multi_class='crammer_singer', random_state=0)
    cs_clf.fit(iris.data, iris.target)

    # similar prediction for ovr and crammer-singer:
    assert_true((ovr_clf.predict(iris.data) ==
                 cs_clf.predict(iris.data)).mean() > .9)

    # classifiers shouldn't be the same
    assert_true((ovr_clf.coef_ != cs_clf.coef_).all())

    # test decision function
    assert_array_equal(cs_clf.predict(iris.data),
                       np.argmax(cs_clf.decision_function(iris.data), axis=1))
    dec_func = np.dot(iris.data, cs_clf.coef_.T) + cs_clf.intercept_
    assert_array_almost_equal(dec_func, cs_clf.decision_function(iris.data))


def test_crammer_singer_binary():
    # Test Crammer-Singer formulation in the binary case
    X, y = make_classification(n_classes=2, random_state=0)

    for fit_intercept in (True, False):
        acc = svm.LinearSVC(fit_intercept=fit_intercept,
                            multi_class="crammer_singer",
                            random_state=0).fit(X, y).score(X, y)
        assert_greater(acc, 0.9)


def test_linearsvc_iris():

    # Test that LinearSVC gives plausible predictions on the iris dataset
    # Also, test symbolic class names (classes_).
    target = iris.target_names[iris.target]
    clf = svm.LinearSVC(random_state=0).fit(iris.data, target)
    assert_equal(set(clf.classes_), set(iris.target_names))
    assert_greater(np.mean(clf.predict(iris.data) == target), 0.8)

    dec = clf.decision_function(iris.data)
    pred = iris.target_names[np.argmax(dec, 1)]
    assert_array_equal(pred, clf.predict(iris.data))


def test_dense_liblinear_intercept_handling(classifier=svm.LinearSVC):
    # Test that dense liblinear honours intercept_scaling param
    X = [[2, 1],
         [3, 1],
         [1, 3],
         [2, 3]]
    y = [0, 0, 1, 1]
    clf = classifier(fit_intercept=True, penalty='l1', loss='squared_hinge',
                     dual=False, C=4, tol=1e-7, random_state=0)
    assert_true(clf.intercept_scaling == 1, clf.intercept_scaling)
    assert_true(clf.fit_intercept)

    # when intercept_scaling is low the intercept value is highly "penalized"
    # by regularization
    clf.intercept_scaling = 1
    clf.fit(X, y)
    assert_almost_equal(clf.intercept_, 0, decimal=5)

    # when intercept_scaling is sufficiently high, the intercept value
    # is not affected by regularization
    clf.intercept_scaling = 100
    clf.fit(X, y)
    intercept1 = clf.intercept_
    assert_less(intercept1, -1)

    # when intercept_scaling is sufficiently high, the intercept value
    # doesn't depend on intercept_scaling value
    clf.intercept_scaling = 1000
    clf.fit(X, y)
    intercept2 = clf.intercept_
    assert_array_almost_equal(intercept1, intercept2, decimal=2)


def test_liblinear_set_coef():
    # multi-class case
    clf = svm.LinearSVC().fit(iris.data, iris.target)
    values = clf.decision_function(iris.data)
    clf.coef_ = clf.coef_.copy()
    clf.intercept_ = clf.intercept_.copy()
    values2 = clf.decision_function(iris.data)
    assert_array_almost_equal(values, values2)

    # binary-class case
    X = [[2, 1],
         [3, 1],
         [1, 3],
         [2, 3]]
    y = [0, 0, 1, 1]

    clf = svm.LinearSVC().fit(X, y)
    values = clf.decision_function(X)
    clf.coef_ = clf.coef_.copy()
    clf.intercept_ = clf.intercept_.copy()
    values2 = clf.decision_function(X)
    assert_array_equal(values, values2)


def test_immutable_coef_property():
    # Check that primal coef modification are not silently ignored
    svms = [
        svm.SVC(kernel='linear').fit(iris.data, iris.target),
        svm.NuSVC(kernel='linear').fit(iris.data, iris.target),
        svm.SVR(kernel='linear').fit(iris.data, iris.target),
        svm.NuSVR(kernel='linear').fit(iris.data, iris.target),
        svm.OneClassSVM(kernel='linear').fit(iris.data),
    ]
    for clf in svms:
        assert_raises(AttributeError, clf.__setattr__, 'coef_', np.arange(3))
        assert_raises((RuntimeError, ValueError),
                      clf.coef_.__setitem__, (0, 0), 0)


def test_linearsvc_verbose():
    # stdout: redirect
    import os
    stdout = os.dup(1)  # save original stdout
    os.dup2(os.pipe()[1], 1)  # replace it

    # actual call
    clf = svm.LinearSVC(verbose=1)
    clf.fit(X, Y)

    # stdout: restore
    os.dup2(stdout, 1)  # restore original stdout


def test_svc_clone_with_callable_kernel():
    # create SVM with callable linear kernel, check that results are the same
    # as with built-in linear kernel
    svm_callable = svm.SVC(kernel=lambda x, y: np.dot(x, y.T),
                           probability=True, random_state=0,
                           decision_function_shape='ovr')
    # clone for checking clonability with lambda functions..
    svm_cloned = base.clone(svm_callable)
    svm_cloned.fit(iris.data, iris.target)

    svm_builtin = svm.SVC(kernel='linear', probability=True, random_state=0,
                          decision_function_shape='ovr')
    svm_builtin.fit(iris.data, iris.target)

    assert_array_almost_equal(svm_cloned.dual_coef_,
                              svm_builtin.dual_coef_)
    assert_array_almost_equal(svm_cloned.intercept_,
                              svm_builtin.intercept_)
    assert_array_equal(svm_cloned.predict(iris.data),
                       svm_builtin.predict(iris.data))

    assert_array_almost_equal(svm_cloned.predict_proba(iris.data),
                              svm_builtin.predict_proba(iris.data),
                              decimal=4)
    assert_array_almost_equal(svm_cloned.decision_function(iris.data),
                              svm_builtin.decision_function(iris.data))


def test_svc_bad_kernel():
    svc = svm.SVC(kernel=lambda x, y: x)
    assert_raises(ValueError, svc.fit, X, Y)


def test_timeout():
    a = svm.SVC(kernel=lambda x, y: np.dot(x, y.T), probability=True,
                random_state=0, max_iter=1)
    assert_warns(ConvergenceWarning, a.fit, X, Y)


def test_unfitted():
    X = "foo!"      # input validation not required when SVM not fitted

    clf = svm.SVC()
    assert_raises_regexp(Exception, r".*\bSVC\b.*\bnot\b.*\bfitted\b",
                         clf.predict, X)

    clf = svm.NuSVR()
    assert_raises_regexp(Exception, r".*\bNuSVR\b.*\bnot\b.*\bfitted\b",
                         clf.predict, X)


# ignore convergence warnings from max_iter=1
@ignore_warnings
def test_consistent_proba():
    a = svm.SVC(probability=True, max_iter=1, random_state=0)
    proba_1 = a.fit(X, Y).predict_proba(X)
    a = svm.SVC(probability=True, max_iter=1, random_state=0)
    proba_2 = a.fit(X, Y).predict_proba(X)
    assert_array_almost_equal(proba_1, proba_2)


def test_linear_svc_convergence_warnings():
    # Test that warnings are raised if model does not converge

    lsvc = svm.LinearSVC(max_iter=2, verbose=1)
    assert_warns(ConvergenceWarning, lsvc.fit, X, Y)
    assert_equal(lsvc.n_iter_, 2)


def test_svr_coef_sign():
    # Test that SVR(kernel="linear") has coef_ with the right sign.
    # Non-regression test for #2933.
    X = np.random.RandomState(21).randn(10, 3)
    y = np.random.RandomState(12).randn(10)

    for svr in [svm.SVR(kernel='linear'), svm.NuSVR(kernel='linear'),
                svm.LinearSVR()]:
        svr.fit(X, y)
        assert_array_almost_equal(svr.predict(X),
                                  np.dot(X, svr.coef_.ravel()) + svr.intercept_)


def test_linear_svc_intercept_scaling():
    # Test that the right error message is thrown when intercept_scaling <= 0

    for i in [-1, 0]:
        lsvc = svm.LinearSVC(intercept_scaling=i)
        msg = ('Intercept scaling is %r but needs to be greater than 0.'
               ' To disable fitting an intercept,'
               ' set fit_intercept=False.' % lsvc.intercept_scaling)
        assert_raise_message(ValueError, msg, lsvc.fit, X, Y)


def test_lsvc_intercept_scaling_zero():
    # Test that intercept_scaling is ignored when fit_intercept is False

    lsvc = svm.LinearSVC(fit_intercept=False)
    lsvc.fit(X, Y)
    assert_equal(lsvc.intercept_, 0.)


def test_hasattr_predict_proba():
    # Method must be (un)available before or after fit, switched by
    # `probability` param

    G = svm.SVC(probability=True)
    assert_true(hasattr(G, 'predict_proba'))
    G.fit(iris.data, iris.target)
    assert_true(hasattr(G, 'predict_proba'))

    G = svm.SVC(probability=False)
    assert_false(hasattr(G, 'predict_proba'))
    G.fit(iris.data, iris.target)
    assert_false(hasattr(G, 'predict_proba'))

    # Switching to `probability=True` after fitting should make
    # predict_proba available, but calling it must not work:
    G.probability = True
    assert_true(hasattr(G, 'predict_proba'))
    msg = "predict_proba is not available when fitted with probability=False"
    assert_raise_message(NotFittedError, msg, G.predict_proba, iris.data)


def test_decision_function_shape_two_class():
    for n_classes in [2, 3]:
        X, y = make_blobs(centers=n_classes, random_state=0)
        for estimator in [svm.SVC, svm.NuSVC]:
            clf = OneVsRestClassifier(estimator(
                decision_function_shape="ovr")).fit(X, y)
            assert_equal(len(clf.predict(X)), len(y))
