import numpy as np
import scipy.sparse as sp
from scipy import linalg, optimize

from sklearn.utils.testing import assert_array_equal
from sklearn.utils.testing import assert_array_almost_equal
from sklearn.utils.testing import assert_equal
from sklearn.utils.testing import assert_greater
from sklearn.utils.testing import assert_raises
from sklearn.utils.testing import assert_true
from sklearn.utils.testing import raises

from sklearn.linear_model import logistic
from sklearn import datasets

X = [[-1, 0], [0, 1], [1, 1]]
X_sp = sp.csr_matrix(X)
Y1 = [0, 1, 1]
Y2 = [2, 1, 0]
iris = datasets.load_iris()


def check_predictions(clf, X, y):
    """Check that the model is able to fit the classification data"""
    n_samples = len(y)
    classes = np.unique(y)
    n_classes = classes.shape[0]

    predicted = clf.fit(X, y).predict(X)
    assert_array_equal(clf.classes_, classes)

    assert_equal(predicted.shape, (n_samples,))
    assert_array_equal(predicted, y)

    probabilities = clf.predict_proba(X)
    assert_equal(probabilities.shape, (n_samples, n_classes))
    assert_array_almost_equal(probabilities.sum(axis=1), np.ones(n_samples))
    assert_array_equal(probabilities.argmax(axis=1), y)


def test_predict_2_classes():
    """Simple sanity check on a 2 classes dataset

    Make sure it predicts the correct result on simple datasets.
    """
    check_predictions(logistic.LogisticRegression(random_state=0), X, Y1)
    check_predictions(logistic.LogisticRegression(random_state=0), X_sp, Y1)

    check_predictions(logistic.LogisticRegression(C=100, random_state=0),
                      X, Y1)
    check_predictions(logistic.LogisticRegression(C=100, random_state=0),
                      X_sp, Y1)

    check_predictions(logistic.LogisticRegression(fit_intercept=False,
                                                  random_state=0), X, Y1)
    check_predictions(logistic.LogisticRegression(fit_intercept=False,
                                                  random_state=0), X_sp, Y1)


def test_error():
    """Test for appropriate exception on errors"""
    assert_raises(ValueError, logistic.LogisticRegression(C=-1).fit, X, Y1)


def test_predict_3_classes():
    check_predictions(logistic.LogisticRegression(C=10), X, Y2)
    check_predictions(logistic.LogisticRegression(C=10), X_sp, Y2)


def test_predict_iris():
    """Test logistic regression with the iris dataset"""
    n_samples, n_features = iris.data.shape

    target = iris.target_names[iris.target]
    clf = logistic.LogisticRegression(C=len(iris.data)).fit(iris.data, target)
    assert_array_equal(np.unique(target), clf.classes_)

    pred = clf.predict(iris.data)
    assert_greater(np.mean(pred == target), .95)

    probabilities = clf.predict_proba(iris.data)
    assert_array_almost_equal(probabilities.sum(axis=1), np.ones(n_samples))

    pred = iris.target_names[probabilities.argmax(axis=1)]
    assert_greater(np.mean(pred == target), .95)


def test_sparsify():
    """Test sparsify and densify members."""
    n_samples, n_features = iris.data.shape
    target = iris.target_names[iris.target]
    clf = logistic.LogisticRegression(random_state=0).fit(iris.data, target)

    pred_d_d = clf.decision_function(iris.data)

    clf.sparsify()
    assert_true(sp.issparse(clf.coef_))
    pred_s_d = clf.decision_function(iris.data)

    sp_data = sp.coo_matrix(iris.data)
    pred_s_s = clf.decision_function(sp_data)

    clf.densify()
    pred_d_s = clf.decision_function(sp_data)

    assert_array_almost_equal(pred_d_d, pred_s_d)
    assert_array_almost_equal(pred_d_d, pred_s_s)
    assert_array_almost_equal(pred_d_d, pred_d_s)


def test_inconsistent_input():
    """Test that an exception is raised on inconsistent input"""
    rng = np.random.RandomState(0)
    X_ = rng.random_sample((5, 10))
    y_ = np.ones(X_.shape[0])
    y_[0] = 0

    clf = logistic.LogisticRegression(random_state=0)

    # Wrong dimensions for training data
    y_wrong = y_[:-1]
    assert_raises(ValueError, clf.fit, X, y_wrong)

    # Wrong dimensions for test data
    assert_raises(ValueError, clf.fit(X_, y_).predict,
                  rng.random_sample((3, 12)))


def test_write_parameters():
    """Test that we can write to coef_ and intercept_"""
    #rng = np.random.RandomState(0)
    #X = rng.random_sample((5, 10))
    #y = np.ones(X.shape[0])
    clf = logistic.LogisticRegression(random_state=0)
    clf.fit(X, Y1)
    clf.coef_[:] = 0
    clf.intercept_[:] = 0
    assert_array_equal(clf.decision_function(X), 0)


@raises(ValueError)
def test_nan():
    """Test proper NaN handling.

    Regression test for Issue #252: fit used to go into an infinite loop.
    """
    Xnan = np.array(X, dtype=np.float64)
    Xnan[0, 1] = np.nan
    logistic.LogisticRegression(random_state=0).fit(Xnan, Y1)


def test_consistency_path():
    """Test that the path algorithm is consistent"""
    Cs = np.logspace(0, 4, 10)
    # can't test with fit_intercept=True since LIBLINEAR
    # penalizes the intercept
    for method in ('lbfgs', 'newton-cg', 'liblinear'):
        coefs, Cs = logistic.logistic_regression_path(
            X, Y1, Cs=Cs, fit_intercept=False, gtol=1e-16, solver=method)
        for i, C in enumerate(Cs):
            lr = logistic.LogisticRegression(
                C=C,fit_intercept=False, tol=1e-16)
            lr.fit(X, Y1)
            lr_coef = lr.coef_.ravel()
            assert_array_almost_equal(lr_coef, coefs[i], decimal=1)

    # test for fit_intercept=True
    for method in ('lbfgs', 'newton-cg', 'liblinear'):
        Cs = [1e3]
        coefs, Cs = logistic.logistic_regression_path(
            X, Y1, Cs=Cs, fit_intercept=True, gtol=1e-16, solver=method)
        lr = logistic.LogisticRegression(
            C=Cs[0], fit_intercept=True, tol=1e-16)
        lr.fit(X, Y1)
        lr_coef = np.concatenate([lr.coef_.ravel(), lr.intercept_])
        assert_array_almost_equal(lr_coef, coefs[0], decimal=1)


def test_liblinear_random_state():
    X, y = datasets.make_classification(n_samples=20)
    lr1 = logistic.LogisticRegression(random_state=0)
    lr1.fit(X, y)
    lr2 = logistic.LogisticRegression(random_state=0)
    lr2.fit(X, y)
    assert_array_almost_equal(lr1.coef_, lr2.coef_)


def test__logistic_loss_and_grad():
    X_ref, y = datasets.make_classification(n_samples=20)
    n_features = X_ref.shape[1]

    X_sp = X_ref.copy()
    X_sp[X_sp < .1] = 0
    X_sp = sp.csr_matrix(X_sp)
    for X in (X_ref, X_sp):
        w = np.zeros(n_features)

        # First check that our derivation of the grad is correct
        loss, grad = logistic._logistic_loss_and_grad(w, X, y, alpha=1.)
        approx_grad = optimize.approx_fprime(w,
                            lambda w: logistic._logistic_loss_and_grad(w, X, y,
                                                            alpha=1.)[0],
                            1e-3
                            )
        assert_array_almost_equal(grad, approx_grad, decimal=2)

        # Second check that our intercept implementation is good
        w = np.zeros(n_features + 1)
        loss_interp, grad_interp =logistic._logistic_loss_and_grad_intercept(w,
                                                    X, y, alpha=1.)
        assert_array_almost_equal(loss, loss_interp)

        approx_grad = optimize.approx_fprime(w,
                            lambda w:
                            logistic._logistic_loss_and_grad_intercept(w, X, y,
                                                            alpha=1.)[0],
                            1e-3
                            )
        assert_array_almost_equal(grad_interp, approx_grad, decimal=2)


def test__logistic_loss_grad_hess():
    n_samples, n_features = 100, 5
    X_ref = np.random.randn(n_samples, n_features)
    y = np.sign(X_ref.dot(5 * np.random.randn(n_features)))
    X_ref -= X_ref.mean()
    X_ref /= X_ref.std()
    X_sp = X_ref.copy()
    X_sp[X_sp < .1] = 0
    X_sp = sp.csr_matrix(X_sp)
    for X in (X_ref, X_sp):
        w = .1 * np.ones(n_features)

        # First check that _logistic_loss_grad_hess is consistent
        # with _logistic_loss_and_grad
        loss, grad = logistic._logistic_loss_and_grad(w, X, y, alpha=1.)
        loss_2, grad_2, hess = logistic._logistic_loss_grad_hess(w, X, y,
                                                    alpha=1.)
        assert_array_almost_equal(grad, grad_2)
        # XXX: we should check a few simple properties of our problem, such
        # as the fact that if X=0, the problem is alpha * ||w||**2, so we
        # know the hessian

        # Now check our hessian along the second direction of the grad
        vector = np.zeros_like(grad)
        vector[1] = 1
        hess_col = hess(vector)

        # Computation of the Hessian is particularly fragile to numerical
        # errors when doing simple finite differences. Here we compute the
        # grad along a path in the direction of the vector and then use a
        # least-square regression to estimate the slope
        e = 1e-3
        d_x = np.linspace(-e, e, 30)
        d_grad = np.array([
                logistic._logistic_loss_and_grad(
                    w + t*vector, X, y, alpha=1.)[1]
                for t in d_x
                ])

        d_grad -= d_grad.mean(axis=0)
        approx_hess_col = linalg.lstsq(d_x[:, np.newaxis], d_grad)[0].ravel()

        assert_array_almost_equal(approx_hess_col, hess_col, decimal=3)

        # Second check that our intercept implementation is good
        w = np.zeros(n_features + 1)
        loss_interp, grad_interp = logistic._logistic_loss_and_grad_intercept(w,
                                                    X, y, alpha=1.)
        loss_interp_2, grad_interp_2, hess = \
                    logistic._logistic_loss_grad_hess_intercept(w,
                                                    X, y, alpha=1.)
        assert_array_almost_equal(loss_interp, loss_interp_2)
        assert_array_almost_equal(grad_interp, grad_interp_2)

def test_logistic_cv():
    # test for LogisticRegressionCV object
    n_samples, n_features = 100, 5
    np.random.seed(0)
    X_ref = np.random.randn(n_samples, n_features)
    y = np.sign(X_ref.dot(5 * np.random.randn(n_features)))
    X_ref -= X_ref.mean()
    X_ref /= X_ref.std()
    lr_cv = logistic.LogisticRegressionCV(Cs=[1.], fit_intercept=False)
    lr_cv.fit(X_ref, y)
    lr = logistic.LogisticRegression(C=1, fit_intercept=False)
    lr.fit(X_ref, y)
    assert_array_almost_equal(lr.coef_, lr_cv.coef_, decimal=2)
