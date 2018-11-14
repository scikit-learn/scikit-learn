import os
import sys
import numpy as np
import scipy.sparse as sp
from scipy import linalg, optimize, sparse

import pytest

from sklearn.base import clone
from sklearn.datasets import load_iris, make_classification
from sklearn.metrics import log_loss
from sklearn.metrics.scorer import get_scorer
from sklearn.model_selection import StratifiedKFold
from sklearn.preprocessing import LabelEncoder
from sklearn.utils import compute_class_weight, _IS_32BIT
from sklearn.utils.testing import assert_almost_equal
from sklearn.utils.testing import assert_allclose
from sklearn.utils.testing import assert_array_almost_equal
from sklearn.utils.testing import assert_array_equal
from sklearn.utils.testing import assert_equal
from sklearn.utils.testing import assert_greater
from sklearn.utils.testing import assert_raise_message
from sklearn.utils.testing import assert_raises
from sklearn.utils.testing import assert_warns
from sklearn.utils.testing import ignore_warnings
from sklearn.utils.testing import assert_warns_message
from sklearn.utils.testing import assert_no_warnings
from sklearn.utils.testing import skip_if_no_parallel

from sklearn.exceptions import ConvergenceWarning
from sklearn.exceptions import ChangedBehaviorWarning
from sklearn.linear_model.logistic import (
    LogisticRegression,
    logistic_regression_path, LogisticRegressionCV,
    _logistic_loss_and_grad, _logistic_grad_hess,
    _multinomial_grad_hess, _logistic_loss,
    _log_reg_scoring_path)

X = [[-1, 0], [0, 1], [1, 1]]
X_sp = sp.csr_matrix(X)
Y1 = [0, 1, 1]
Y2 = [2, 1, 0]
iris = load_iris()


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


@pytest.mark.filterwarnings('ignore: Default solver will be changed')  # 0.22
@pytest.mark.filterwarnings('ignore: Default multi_class will')  # 0.22
def test_predict_2_classes():
    # Simple sanity check on a 2 classes dataset
    # Make sure it predicts the correct result on simple datasets.
    check_predictions(LogisticRegression(random_state=0), X, Y1)
    check_predictions(LogisticRegression(random_state=0), X_sp, Y1)

    check_predictions(LogisticRegression(C=100, random_state=0), X, Y1)
    check_predictions(LogisticRegression(C=100, random_state=0), X_sp, Y1)

    check_predictions(LogisticRegression(fit_intercept=False,
                                         random_state=0), X, Y1)
    check_predictions(LogisticRegression(fit_intercept=False,
                                         random_state=0), X_sp, Y1)


@pytest.mark.filterwarnings('ignore: Default multi_class will')  # 0.22
def test_error():
    # Test for appropriate exception on errors
    msg = "Penalty term must be positive"
    assert_raise_message(ValueError, msg,
                         LogisticRegression(C=-1).fit, X, Y1)
    assert_raise_message(ValueError, msg,
                         LogisticRegression(C="test").fit, X, Y1)

    msg = "is not a valid scoring value"
    assert_raise_message(ValueError, msg,
                         LogisticRegressionCV(scoring='bad-scorer', cv=2).fit,
                         X, Y1)

    for LR in [LogisticRegression, LogisticRegressionCV]:
        msg = "Tolerance for stopping criteria must be positive"
        assert_raise_message(ValueError, msg, LR(tol=-1).fit, X, Y1)
        assert_raise_message(ValueError, msg, LR(tol="test").fit, X, Y1)

        msg = "Maximum number of iteration must be positive"
        assert_raise_message(ValueError, msg, LR(max_iter=-1).fit, X, Y1)
        assert_raise_message(ValueError, msg, LR(max_iter="test").fit, X, Y1)


@pytest.mark.filterwarnings('ignore: Default multi_class will')  # 0.22
def test_logistic_cv_mock_scorer():

    class MockScorer(object):
        def __init__(self):
            self.calls = 0
            self.scores = [0.1, 0.4, 0.8, 0.5]

        def __call__(self, model, X, y, sample_weight=None):
            score = self.scores[self.calls % len(self.scores)]
            self.calls += 1
            return score

    mock_scorer = MockScorer()
    Cs = [1, 2, 3, 4]
    cv = 2

    lr = LogisticRegressionCV(Cs=Cs, scoring=mock_scorer, cv=cv)
    lr.fit(X, Y1)

    # Cs[2] has the highest score (0.8) from MockScorer
    assert lr.C_[0] == Cs[2]

    # scorer called 8 times (cv*len(Cs))
    assert mock_scorer.calls == cv * len(Cs)

    # reset mock_scorer
    mock_scorer.calls = 0
    with pytest.warns(ChangedBehaviorWarning):
        custom_score = lr.score(X, lr.predict(X))

    assert custom_score == mock_scorer.scores[0]
    assert mock_scorer.calls == 1


def test_logistic_cv_score_does_not_warn_by_default():
    lr = LogisticRegressionCV(cv=2, multi_class='ovr')
    lr.fit(X, Y1)

    with pytest.warns(None) as record:
        lr.score(X, lr.predict(X))
    assert len(record) == 0


@skip_if_no_parallel
def test_lr_liblinear_warning():
    n_samples, n_features = iris.data.shape
    target = iris.target_names[iris.target]

    lr = LogisticRegression(solver='liblinear', multi_class='ovr', n_jobs=2)
    assert_warns_message(UserWarning,
                         "'n_jobs' > 1 does not have any effect when"
                         " 'solver' is set to 'liblinear'. Got 'n_jobs'"
                         " = 2.",
                         lr.fit, iris.data, target)


@pytest.mark.filterwarnings('ignore: Default solver will be changed')  # 0.22
@pytest.mark.filterwarnings('ignore: Default multi_class will')  # 0.22
def test_predict_3_classes():
    check_predictions(LogisticRegression(C=10), X, Y2)
    check_predictions(LogisticRegression(C=10), X_sp, Y2)


def test_predict_iris():
    # Test logistic regression with the iris dataset
    n_samples, n_features = iris.data.shape

    target = iris.target_names[iris.target]

    # Test that both multinomial and OvR solvers handle
    # multiclass data correctly and give good accuracy
    # score (>0.95) for the training data.
    for clf in [LogisticRegression(C=len(iris.data), solver='liblinear',
                                   multi_class='ovr'),
                LogisticRegression(C=len(iris.data), solver='lbfgs',
                                   multi_class='multinomial'),
                LogisticRegression(C=len(iris.data), solver='newton-cg',
                                   multi_class='multinomial'),
                LogisticRegression(C=len(iris.data), solver='sag', tol=1e-2,
                                   multi_class='ovr', random_state=42),
                LogisticRegression(C=len(iris.data), solver='saga', tol=1e-2,
                                   multi_class='ovr', random_state=42)
                ]:
        clf.fit(iris.data, target)
        assert_array_equal(np.unique(target), clf.classes_)

        pred = clf.predict(iris.data)
        assert_greater(np.mean(pred == target), .95)

        probabilities = clf.predict_proba(iris.data)
        assert_array_almost_equal(probabilities.sum(axis=1),
                                  np.ones(n_samples))

        pred = iris.target_names[probabilities.argmax(axis=1)]
        assert_greater(np.mean(pred == target), .95)


@pytest.mark.parametrize('solver', ['lbfgs', 'newton-cg', 'sag', 'saga'])
def test_multinomial_validation(solver):
    lr = LogisticRegression(C=-1, solver=solver, multi_class='multinomial')
    assert_raises(ValueError, lr.fit, [[0, 1], [1, 0]], [0, 1])


@pytest.mark.parametrize('LR', [LogisticRegression, LogisticRegressionCV])
def test_check_solver_option(LR):
    X, y = iris.data, iris.target

    msg = ("Logistic Regression supports only solvers in ['liblinear', "
           "'newton-cg', 'lbfgs', 'sag', 'saga'], got wrong_name.")
    lr = LR(solver="wrong_name", multi_class="ovr")
    assert_raise_message(ValueError, msg, lr.fit, X, y)

    msg = ("multi_class should be 'multinomial', 'ovr' or 'auto'. "
           "Got wrong_name")
    lr = LR(solver='newton-cg', multi_class="wrong_name")
    assert_raise_message(ValueError, msg, lr.fit, X, y)

    # only 'liblinear' solver
    msg = "Solver liblinear does not support a multinomial backend."
    lr = LR(solver='liblinear', multi_class='multinomial')
    assert_raise_message(ValueError, msg, lr.fit, X, y)

    # all solvers except 'liblinear'
    for solver in ['newton-cg', 'lbfgs', 'sag']:
        msg = ("Solver %s supports only l2 penalties, got l1 penalty." %
               solver)
        lr = LR(solver=solver, penalty='l1', multi_class='ovr')
        assert_raise_message(ValueError, msg, lr.fit, X, y)
    for solver in ['newton-cg', 'lbfgs', 'sag', 'saga']:
        msg = ("Solver %s supports only dual=False, got dual=True" %
               solver)
        lr = LR(solver=solver, dual=True, multi_class='ovr')
        assert_raise_message(ValueError, msg, lr.fit, X, y)


@pytest.mark.parametrize('model, params, warn_solver',
                         [(LogisticRegression, {}, True),
                          (LogisticRegressionCV, {'cv': 5}, False)])
def test_logistic_regression_warnings(model, params, warn_solver):
    clf_solver_warning = model(multi_class='ovr', **params)
    clf_multi_class_warning = model(solver='lbfgs', **params)
    clf_no_warnings = model(solver='lbfgs', multi_class='ovr', **params)

    solver_warning_msg = "Default solver will be changed to 'lbfgs'"
    multi_class_warning_msg = "Default multi_class will be changed to 'auto"

    if warn_solver:
        assert_warns_message(FutureWarning, solver_warning_msg,
                             clf_solver_warning.fit, iris.data, iris.target)
    else:
        assert_no_warnings(clf_no_warnings.fit, iris.data, iris.target)

    assert_warns_message(FutureWarning, multi_class_warning_msg,
                         clf_multi_class_warning.fit, iris.data, iris.target)
    # But no warning when binary target:
    assert_no_warnings(clf_multi_class_warning.fit,
                       iris.data, iris.target == 0)
    assert_no_warnings(clf_no_warnings.fit, iris.data, iris.target)


@pytest.mark.parametrize('solver', ['lbfgs', 'newton-cg', 'sag', 'saga'])
def test_multinomial_binary(solver):
    # Test multinomial LR on a binary problem.
    target = (iris.target > 0).astype(np.intp)
    target = np.array(["setosa", "not-setosa"])[target]

    clf = LogisticRegression(solver=solver, multi_class='multinomial',
                             random_state=42, max_iter=2000)
    clf.fit(iris.data, target)

    assert_equal(clf.coef_.shape, (1, iris.data.shape[1]))
    assert_equal(clf.intercept_.shape, (1,))
    assert_array_equal(clf.predict(iris.data), target)

    mlr = LogisticRegression(solver=solver, multi_class='multinomial',
                             random_state=42, fit_intercept=False)
    mlr.fit(iris.data, target)
    pred = clf.classes_[np.argmax(clf.predict_log_proba(iris.data),
                                  axis=1)]
    assert_greater(np.mean(pred == target), .9)


def test_multinomial_binary_probabilities():
    # Test multinomial LR gives expected probabilities based on the
    # decision function, for a binary problem.
    X, y = make_classification()
    clf = LogisticRegression(multi_class='multinomial', solver='saga')
    clf.fit(X, y)

    decision = clf.decision_function(X)
    proba = clf.predict_proba(X)

    expected_proba_class_1 = (np.exp(decision) /
                              (np.exp(decision) + np.exp(-decision)))
    expected_proba = np.c_[1 - expected_proba_class_1, expected_proba_class_1]

    assert_almost_equal(proba, expected_proba)


@pytest.mark.filterwarnings('ignore: Default solver will be changed')  # 0.22
@pytest.mark.filterwarnings('ignore: Default multi_class will')  # 0.22
def test_sparsify():
    # Test sparsify and densify members.
    n_samples, n_features = iris.data.shape
    target = iris.target_names[iris.target]
    clf = LogisticRegression(random_state=0).fit(iris.data, target)

    pred_d_d = clf.decision_function(iris.data)

    clf.sparsify()
    assert sp.issparse(clf.coef_)
    pred_s_d = clf.decision_function(iris.data)

    sp_data = sp.coo_matrix(iris.data)
    pred_s_s = clf.decision_function(sp_data)

    clf.densify()
    pred_d_s = clf.decision_function(sp_data)

    assert_array_almost_equal(pred_d_d, pred_s_d)
    assert_array_almost_equal(pred_d_d, pred_s_s)
    assert_array_almost_equal(pred_d_d, pred_d_s)


@pytest.mark.filterwarnings('ignore: Default solver will be changed')  # 0.22
@pytest.mark.filterwarnings('ignore: Default multi_class will')  # 0.22
def test_inconsistent_input():
    # Test that an exception is raised on inconsistent input
    rng = np.random.RandomState(0)
    X_ = rng.random_sample((5, 10))
    y_ = np.ones(X_.shape[0])
    y_[0] = 0

    clf = LogisticRegression(random_state=0)

    # Wrong dimensions for training data
    y_wrong = y_[:-1]
    assert_raises(ValueError, clf.fit, X, y_wrong)

    # Wrong dimensions for test data
    assert_raises(ValueError, clf.fit(X_, y_).predict,
                  rng.random_sample((3, 12)))


@pytest.mark.filterwarnings('ignore: Default solver will be changed')  # 0.22
@pytest.mark.filterwarnings('ignore: Default multi_class will')  # 0.22
def test_write_parameters():
    # Test that we can write to coef_ and intercept_
    clf = LogisticRegression(random_state=0)
    clf.fit(X, Y1)
    clf.coef_[:] = 0
    clf.intercept_[:] = 0
    assert_array_almost_equal(clf.decision_function(X), 0)


@pytest.mark.filterwarnings('ignore: Default solver will be changed')  # 0.22
@pytest.mark.filterwarnings('ignore: Default multi_class will')  # 0.22
def test_nan():
    # Test proper NaN handling.
    # Regression test for Issue #252: fit used to go into an infinite loop.
    Xnan = np.array(X, dtype=np.float64)
    Xnan[0, 1] = np.nan
    logistic = LogisticRegression(random_state=0)
    assert_raises(ValueError, logistic.fit, Xnan, Y1)


def test_consistency_path():
    # Test that the path algorithm is consistent
    rng = np.random.RandomState(0)
    X = np.concatenate((rng.randn(100, 2) + [1, 1], rng.randn(100, 2)))
    y = [1] * 100 + [-1] * 100
    Cs = np.logspace(0, 4, 10)

    f = ignore_warnings
    # can't test with fit_intercept=True since LIBLINEAR
    # penalizes the intercept
    for solver in ['sag', 'saga']:
        coefs, Cs, _ = f(logistic_regression_path)(
            X, y, Cs=Cs, fit_intercept=False, tol=1e-5, solver=solver,
            max_iter=1000, multi_class='ovr', random_state=0)
        for i, C in enumerate(Cs):
            lr = LogisticRegression(C=C, fit_intercept=False, tol=1e-5,
                                    solver=solver, multi_class='ovr',
                                    random_state=0, max_iter=1000)
            lr.fit(X, y)
            lr_coef = lr.coef_.ravel()
            assert_array_almost_equal(lr_coef, coefs[i], decimal=4,
                                      err_msg="with solver = %s" % solver)

    # test for fit_intercept=True
    for solver in ('lbfgs', 'newton-cg', 'liblinear', 'sag', 'saga'):
        Cs = [1e3]
        coefs, Cs, _ = f(logistic_regression_path)(
            X, y, Cs=Cs, fit_intercept=True, tol=1e-6, solver=solver,
            intercept_scaling=10000., random_state=0, multi_class='ovr')
        lr = LogisticRegression(C=Cs[0], fit_intercept=True, tol=1e-4,
                                intercept_scaling=10000., random_state=0,
                                multi_class='ovr', solver=solver)
        lr.fit(X, y)
        lr_coef = np.concatenate([lr.coef_.ravel(), lr.intercept_])
        assert_array_almost_equal(lr_coef, coefs[0], decimal=4,
                                  err_msg="with solver = %s" % solver)


def test_logistic_regression_path_convergence_fail():
    rng = np.random.RandomState(0)
    X = np.concatenate((rng.randn(100, 2) + [1, 1], rng.randn(100, 2)))
    y = [1] * 100 + [-1] * 100
    Cs = [1e3]
    assert_warns(ConvergenceWarning, logistic_regression_path,
                 X, y, Cs=Cs, tol=0., max_iter=1, random_state=0, verbose=1)


def test_liblinear_dual_random_state():
    # random_state is relevant for liblinear solver only if dual=True
    X, y = make_classification(n_samples=20, random_state=0)
    lr1 = LogisticRegression(random_state=0, dual=True, max_iter=1, tol=1e-15,
                             solver='liblinear', multi_class='ovr')
    lr1.fit(X, y)
    lr2 = LogisticRegression(random_state=0, dual=True, max_iter=1, tol=1e-15,
                             solver='liblinear', multi_class='ovr')
    lr2.fit(X, y)
    lr3 = LogisticRegression(random_state=8, dual=True, max_iter=1, tol=1e-15,
                             solver='liblinear', multi_class='ovr')
    lr3.fit(X, y)

    # same result for same random state
    assert_array_almost_equal(lr1.coef_, lr2.coef_)
    # different results for different random states
    msg = "Arrays are not almost equal to 6 decimals"
    assert_raise_message(AssertionError, msg,
                         assert_array_almost_equal, lr1.coef_, lr3.coef_)


def test_logistic_loss_and_grad():
    X_ref, y = make_classification(n_samples=20, random_state=0)
    n_features = X_ref.shape[1]

    X_sp = X_ref.copy()
    X_sp[X_sp < .1] = 0
    X_sp = sp.csr_matrix(X_sp)
    for X in (X_ref, X_sp):
        w = np.zeros(n_features)

        # First check that our derivation of the grad is correct
        loss, grad = _logistic_loss_and_grad(w, X, y, alpha=1.)
        approx_grad = optimize.approx_fprime(
            w, lambda w: _logistic_loss_and_grad(w, X, y, alpha=1.)[0], 1e-3
        )
        assert_array_almost_equal(grad, approx_grad, decimal=2)

        # Second check that our intercept implementation is good
        w = np.zeros(n_features + 1)
        loss_interp, grad_interp = _logistic_loss_and_grad(
            w, X, y, alpha=1.
        )
        assert_array_almost_equal(loss, loss_interp)

        approx_grad = optimize.approx_fprime(
            w, lambda w: _logistic_loss_and_grad(w, X, y, alpha=1.)[0], 1e-3
        )
        assert_array_almost_equal(grad_interp, approx_grad, decimal=2)


def test_logistic_grad_hess():
    rng = np.random.RandomState(0)
    n_samples, n_features = 50, 5
    X_ref = rng.randn(n_samples, n_features)
    y = np.sign(X_ref.dot(5 * rng.randn(n_features)))
    X_ref -= X_ref.mean()
    X_ref /= X_ref.std()
    X_sp = X_ref.copy()
    X_sp[X_sp < .1] = 0
    X_sp = sp.csr_matrix(X_sp)
    for X in (X_ref, X_sp):
        w = np.full(n_features, .1)

        # First check that _logistic_grad_hess is consistent
        # with _logistic_loss_and_grad
        loss, grad = _logistic_loss_and_grad(w, X, y, alpha=1.)
        grad_2, hess = _logistic_grad_hess(w, X, y, alpha=1.)
        assert_array_almost_equal(grad, grad_2)

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
            _logistic_loss_and_grad(w + t * vector, X, y, alpha=1.)[1]
            for t in d_x
        ])

        d_grad -= d_grad.mean(axis=0)
        approx_hess_col = linalg.lstsq(d_x[:, np.newaxis], d_grad)[0].ravel()

        assert_array_almost_equal(approx_hess_col, hess_col, decimal=3)

        # Second check that our intercept implementation is good
        w = np.zeros(n_features + 1)
        loss_interp, grad_interp = _logistic_loss_and_grad(w, X, y, alpha=1.)
        loss_interp_2 = _logistic_loss(w, X, y, alpha=1.)
        grad_interp_2, hess = _logistic_grad_hess(w, X, y, alpha=1.)
        assert_array_almost_equal(loss_interp, loss_interp_2)
        assert_array_almost_equal(grad_interp, grad_interp_2)


@pytest.mark.filterwarnings('ignore: You should specify a value')  # 0.22
def test_logistic_cv():
    # test for LogisticRegressionCV object
    n_samples, n_features = 50, 5
    rng = np.random.RandomState(0)
    X_ref = rng.randn(n_samples, n_features)
    y = np.sign(X_ref.dot(5 * rng.randn(n_features)))
    X_ref -= X_ref.mean()
    X_ref /= X_ref.std()
    lr_cv = LogisticRegressionCV(Cs=[1.], fit_intercept=False,
                                 solver='liblinear', multi_class='ovr')
    lr_cv.fit(X_ref, y)
    lr = LogisticRegression(C=1., fit_intercept=False,
                            solver='liblinear', multi_class='ovr')
    lr.fit(X_ref, y)
    assert_array_almost_equal(lr.coef_, lr_cv.coef_)

    assert_array_equal(lr_cv.coef_.shape, (1, n_features))
    assert_array_equal(lr_cv.classes_, [-1, 1])
    assert_equal(len(lr_cv.classes_), 2)

    coefs_paths = np.asarray(list(lr_cv.coefs_paths_.values()))
    assert_array_equal(coefs_paths.shape, (1, 3, 1, n_features))
    assert_array_equal(lr_cv.Cs_.shape, (1,))
    scores = np.asarray(list(lr_cv.scores_.values()))
    assert_array_equal(scores.shape, (1, 3, 1))


@pytest.mark.filterwarnings('ignore: You should specify a value')  # 0.22
@pytest.mark.parametrize('scoring, multiclass_agg_list',
                         [('accuracy', ['']),
                          ('precision', ['_macro', '_weighted']),
                          # no need to test for micro averaging because it
                          # is the same as accuracy for f1, precision,
                          # and recall (see https://github.com/
                          # scikit-learn/scikit-learn/pull/
                          # 11578#discussion_r203250062)
                          ('f1', ['_macro', '_weighted']),
                          ('neg_log_loss', ['']),
                          ('recall', ['_macro', '_weighted'])])
def test_logistic_cv_multinomial_score(scoring, multiclass_agg_list):
    # test that LogisticRegressionCV uses the right score to compute its
    # cross-validation scores when using a multinomial scoring
    # see https://github.com/scikit-learn/scikit-learn/issues/8720
    X, y = make_classification(n_samples=100, random_state=0, n_classes=3,
                               n_informative=6)
    train, test = np.arange(80), np.arange(80, 100)
    lr = LogisticRegression(C=1., solver='lbfgs', multi_class='multinomial')
    # we use lbfgs to support multinomial
    params = lr.get_params()
    # we store the params to set them further in _log_reg_scoring_path
    for key in ['C', 'n_jobs', 'warm_start']:
        del params[key]
    lr.fit(X[train], y[train])
    for averaging in multiclass_agg_list:
        scorer = get_scorer(scoring + averaging)
        assert_array_almost_equal(
            _log_reg_scoring_path(X, y, train, test, Cs=[1.],
                                  scoring=scorer, **params)[2][0],
            scorer(lr, X[test], y[test]))


@pytest.mark.filterwarnings('ignore: You should specify a value')  # 0.22
def test_multinomial_logistic_regression_string_inputs():
    # Test with string labels for LogisticRegression(CV)
    n_samples, n_features, n_classes = 50, 5, 3
    X_ref, y = make_classification(n_samples=n_samples, n_features=n_features,
                                   n_classes=n_classes, n_informative=3,
                                   random_state=0)
    y_str = LabelEncoder().fit(['bar', 'baz', 'foo']).inverse_transform(y)
    # For numerical labels, let y values be taken from set (-1, 0, 1)
    y = np.array(y) - 1
    # Test for string labels
    lr = LogisticRegression(solver='lbfgs', multi_class='multinomial')
    lr_cv = LogisticRegressionCV(solver='lbfgs', multi_class='multinomial')
    lr_str = LogisticRegression(solver='lbfgs', multi_class='multinomial')
    lr_cv_str = LogisticRegressionCV(solver='lbfgs', multi_class='multinomial')

    lr.fit(X_ref, y)
    lr_cv.fit(X_ref, y)
    lr_str.fit(X_ref, y_str)
    lr_cv_str.fit(X_ref, y_str)

    assert_array_almost_equal(lr.coef_, lr_str.coef_)
    assert_equal(sorted(lr_str.classes_), ['bar', 'baz', 'foo'])
    assert_array_almost_equal(lr_cv.coef_, lr_cv_str.coef_)
    assert_equal(sorted(lr_str.classes_), ['bar', 'baz', 'foo'])
    assert_equal(sorted(lr_cv_str.classes_), ['bar', 'baz', 'foo'])

    # The predictions should be in original labels
    assert_equal(sorted(np.unique(lr_str.predict(X_ref))),
                 ['bar', 'baz', 'foo'])
    assert_equal(sorted(np.unique(lr_cv_str.predict(X_ref))),
                 ['bar', 'baz', 'foo'])

    # Make sure class weights can be given with string labels
    lr_cv_str = LogisticRegression(
        solver='lbfgs', class_weight={'bar': 1, 'baz': 2, 'foo': 0},
        multi_class='multinomial').fit(X_ref, y_str)
    assert_equal(sorted(np.unique(lr_cv_str.predict(X_ref))), ['bar', 'baz'])


@pytest.mark.filterwarnings('ignore: Default multi_class will')  # 0.22
@pytest.mark.filterwarnings('ignore: You should specify a value')  # 0.22
def test_logistic_cv_sparse():
    X, y = make_classification(n_samples=50, n_features=5,
                               random_state=0)
    X[X < 1.0] = 0.0
    csr = sp.csr_matrix(X)

    clf = LogisticRegressionCV(fit_intercept=True)
    clf.fit(X, y)
    clfs = LogisticRegressionCV(fit_intercept=True)
    clfs.fit(csr, y)
    assert_array_almost_equal(clfs.coef_, clf.coef_)
    assert_array_almost_equal(clfs.intercept_, clf.intercept_)
    assert_equal(clfs.C_, clf.C_)


def test_intercept_logistic_helper():
    n_samples, n_features = 10, 5
    X, y = make_classification(n_samples=n_samples, n_features=n_features,
                               random_state=0)

    # Fit intercept case.
    alpha = 1.
    w = np.ones(n_features + 1)
    grad_interp, hess_interp = _logistic_grad_hess(w, X, y, alpha)
    loss_interp = _logistic_loss(w, X, y, alpha)

    # Do not fit intercept. This can be considered equivalent to adding
    # a feature vector of ones, i.e column of one vectors.
    X_ = np.hstack((X, np.ones(10)[:, np.newaxis]))
    grad, hess = _logistic_grad_hess(w, X_, y, alpha)
    loss = _logistic_loss(w, X_, y, alpha)

    # In the fit_intercept=False case, the feature vector of ones is
    # penalized. This should be taken care of.
    assert_almost_equal(loss_interp + 0.5 * (w[-1] ** 2), loss)

    # Check gradient.
    assert_array_almost_equal(grad_interp[:n_features], grad[:n_features])
    assert_almost_equal(grad_interp[-1] + alpha * w[-1], grad[-1])

    rng = np.random.RandomState(0)
    grad = rng.rand(n_features + 1)
    hess_interp = hess_interp(grad)
    hess = hess(grad)
    assert_array_almost_equal(hess_interp[:n_features], hess[:n_features])
    assert_almost_equal(hess_interp[-1] + alpha * grad[-1], hess[-1])


def test_ovr_multinomial_iris():
    # Test that OvR and multinomial are correct using the iris dataset.
    train, target = iris.data, iris.target
    n_samples, n_features = train.shape

    # The cv indices from stratified kfold (where stratification is done based
    # on the fine-grained iris classes, i.e, before the classes 0 and 1 are
    # conflated) is used for both clf and clf1
    n_cv = 2
    cv = StratifiedKFold(n_cv)
    precomputed_folds = list(cv.split(train, target))

    # Train clf on the original dataset where classes 0 and 1 are separated
    clf = LogisticRegressionCV(cv=precomputed_folds, multi_class='ovr')
    clf.fit(train, target)

    # Conflate classes 0 and 1 and train clf1 on this modified dataset
    clf1 = LogisticRegressionCV(cv=precomputed_folds, multi_class='ovr')
    target_copy = target.copy()
    target_copy[target_copy == 0] = 1
    clf1.fit(train, target_copy)

    # Ensure that what OvR learns for class2 is same regardless of whether
    # classes 0 and 1 are separated or not
    assert_array_almost_equal(clf.scores_[2], clf1.scores_[2])
    assert_array_almost_equal(clf.intercept_[2:], clf1.intercept_)
    assert_array_almost_equal(clf.coef_[2][np.newaxis, :], clf1.coef_)

    # Test the shape of various attributes.
    assert_equal(clf.coef_.shape, (3, n_features))
    assert_array_equal(clf.classes_, [0, 1, 2])
    coefs_paths = np.asarray(list(clf.coefs_paths_.values()))
    assert_array_almost_equal(coefs_paths.shape, (3, n_cv, 10, n_features + 1))
    assert_equal(clf.Cs_.shape, (10,))
    scores = np.asarray(list(clf.scores_.values()))
    assert_equal(scores.shape, (3, n_cv, 10))

    # Test that for the iris data multinomial gives a better accuracy than OvR
    for solver in ['lbfgs', 'newton-cg', 'sag', 'saga']:
        max_iter = 2000 if solver in ['sag', 'saga'] else 15
        clf_multi = LogisticRegressionCV(
            solver=solver, multi_class='multinomial', max_iter=max_iter,
            random_state=42, tol=1e-5 if solver in ['sag', 'saga'] else 1e-2,
            cv=2)
        clf_multi.fit(train, target)
        multi_score = clf_multi.score(train, target)
        ovr_score = clf.score(train, target)
        assert_greater(multi_score, ovr_score)

        # Test attributes of LogisticRegressionCV
        assert_equal(clf.coef_.shape, clf_multi.coef_.shape)
        assert_array_equal(clf_multi.classes_, [0, 1, 2])
        coefs_paths = np.asarray(list(clf_multi.coefs_paths_.values()))
        assert_array_almost_equal(coefs_paths.shape, (3, n_cv, 10,
                                                      n_features + 1))
        assert_equal(clf_multi.Cs_.shape, (10,))
        scores = np.asarray(list(clf_multi.scores_.values()))
        assert_equal(scores.shape, (3, n_cv, 10))


def test_logistic_regression_solvers():
    X, y = make_classification(n_features=10, n_informative=5, random_state=0)

    params = dict(fit_intercept=False, random_state=42, multi_class='ovr')
    ncg = LogisticRegression(solver='newton-cg', **params)
    lbf = LogisticRegression(solver='lbfgs', **params)
    lib = LogisticRegression(solver='liblinear', **params)
    sag = LogisticRegression(solver='sag', **params)
    saga = LogisticRegression(solver='saga', **params)
    ncg.fit(X, y)
    lbf.fit(X, y)
    sag.fit(X, y)
    saga.fit(X, y)
    lib.fit(X, y)
    assert_array_almost_equal(ncg.coef_, lib.coef_, decimal=3)
    assert_array_almost_equal(lib.coef_, lbf.coef_, decimal=3)
    assert_array_almost_equal(ncg.coef_, lbf.coef_, decimal=3)
    assert_array_almost_equal(sag.coef_, lib.coef_, decimal=3)
    assert_array_almost_equal(sag.coef_, ncg.coef_, decimal=3)
    assert_array_almost_equal(sag.coef_, lbf.coef_, decimal=3)
    assert_array_almost_equal(saga.coef_, sag.coef_, decimal=3)
    assert_array_almost_equal(saga.coef_, lbf.coef_, decimal=3)
    assert_array_almost_equal(saga.coef_, ncg.coef_, decimal=3)
    assert_array_almost_equal(saga.coef_, lib.coef_, decimal=3)


def test_logistic_regression_solvers_multiclass():
    X, y = make_classification(n_samples=20, n_features=20, n_informative=10,
                               n_classes=3, random_state=0)
    tol = 1e-7
    params = dict(fit_intercept=False, tol=tol, random_state=42,
                  multi_class='ovr')
    ncg = LogisticRegression(solver='newton-cg', **params)
    lbf = LogisticRegression(solver='lbfgs', **params)
    lib = LogisticRegression(solver='liblinear', **params)
    sag = LogisticRegression(solver='sag', max_iter=1000, **params)
    saga = LogisticRegression(solver='saga', max_iter=10000, **params)
    ncg.fit(X, y)
    lbf.fit(X, y)
    sag.fit(X, y)
    saga.fit(X, y)
    lib.fit(X, y)
    assert_array_almost_equal(ncg.coef_, lib.coef_, decimal=4)
    assert_array_almost_equal(lib.coef_, lbf.coef_, decimal=4)
    assert_array_almost_equal(ncg.coef_, lbf.coef_, decimal=4)
    assert_array_almost_equal(sag.coef_, lib.coef_, decimal=4)
    assert_array_almost_equal(sag.coef_, ncg.coef_, decimal=4)
    assert_array_almost_equal(sag.coef_, lbf.coef_, decimal=4)
    assert_array_almost_equal(saga.coef_, sag.coef_, decimal=4)
    assert_array_almost_equal(saga.coef_, lbf.coef_, decimal=4)
    assert_array_almost_equal(saga.coef_, ncg.coef_, decimal=4)
    assert_array_almost_equal(saga.coef_, lib.coef_, decimal=4)


@pytest.mark.filterwarnings('ignore: You should specify a value')  # 0.22
def test_logistic_regressioncv_class_weights():
    for weight in [{0: 0.1, 1: 0.2}, {0: 0.1, 1: 0.2, 2: 0.5}]:
        n_classes = len(weight)
        for class_weight in (weight, 'balanced'):
            X, y = make_classification(n_samples=30, n_features=3,
                                       n_repeated=0,
                                       n_informative=3, n_redundant=0,
                                       n_classes=n_classes, random_state=0)

            clf_lbf = LogisticRegressionCV(solver='lbfgs', Cs=1,
                                           fit_intercept=False,
                                           multi_class='ovr',
                                           class_weight=class_weight)
            clf_ncg = LogisticRegressionCV(solver='newton-cg', Cs=1,
                                           fit_intercept=False,
                                           multi_class='ovr',
                                           class_weight=class_weight)
            clf_lib = LogisticRegressionCV(solver='liblinear', Cs=1,
                                           fit_intercept=False,
                                           multi_class='ovr',
                                           class_weight=class_weight)
            clf_sag = LogisticRegressionCV(solver='sag', Cs=1,
                                           fit_intercept=False,
                                           multi_class='ovr',
                                           class_weight=class_weight,
                                           tol=1e-5, max_iter=10000,
                                           random_state=0)
            clf_saga = LogisticRegressionCV(solver='saga', Cs=1,
                                            fit_intercept=False,
                                            multi_class='ovr',
                                            class_weight=class_weight,
                                            tol=1e-5, max_iter=10000,
                                            random_state=0)
            clf_lbf.fit(X, y)
            clf_ncg.fit(X, y)
            clf_lib.fit(X, y)
            clf_sag.fit(X, y)
            clf_saga.fit(X, y)
            assert_array_almost_equal(clf_lib.coef_, clf_lbf.coef_, decimal=4)
            assert_array_almost_equal(clf_ncg.coef_, clf_lbf.coef_, decimal=4)
            assert_array_almost_equal(clf_sag.coef_, clf_lbf.coef_, decimal=4)
            assert_array_almost_equal(clf_saga.coef_, clf_lbf.coef_, decimal=4)


@pytest.mark.filterwarnings('ignore: You should specify a value')  # 0.22
def test_logistic_regression_sample_weights():
    X, y = make_classification(n_samples=20, n_features=5, n_informative=3,
                               n_classes=2, random_state=0)
    sample_weight = y + 1

    for LR in [LogisticRegression, LogisticRegressionCV]:

        # Test that passing sample_weight as ones is the same as
        # not passing them at all (default None)
        for solver in ['lbfgs', 'liblinear']:
            clf_sw_none = LR(solver=solver, fit_intercept=False,
                             random_state=42, multi_class='ovr')
            clf_sw_none.fit(X, y)
            clf_sw_ones = LR(solver=solver, fit_intercept=False,
                             random_state=42, multi_class='ovr')
            clf_sw_ones.fit(X, y, sample_weight=np.ones(y.shape[0]))
            assert_array_almost_equal(
                clf_sw_none.coef_, clf_sw_ones.coef_, decimal=4)

        # Test that sample weights work the same with the lbfgs,
        # newton-cg, and 'sag' solvers
        clf_sw_lbfgs = LR(solver='lbfgs', fit_intercept=False, random_state=42,
                          multi_class='ovr')
        clf_sw_lbfgs.fit(X, y, sample_weight=sample_weight)
        clf_sw_n = LR(solver='newton-cg', fit_intercept=False, random_state=42,
                      multi_class='ovr')
        clf_sw_n.fit(X, y, sample_weight=sample_weight)
        clf_sw_sag = LR(solver='sag', fit_intercept=False, tol=1e-10,
                        random_state=42, multi_class='ovr')
        # ignore convergence warning due to small dataset
        with ignore_warnings():
            clf_sw_sag.fit(X, y, sample_weight=sample_weight)
        clf_sw_liblinear = LR(solver='liblinear', fit_intercept=False,
                              random_state=42, multi_class='ovr')
        clf_sw_liblinear.fit(X, y, sample_weight=sample_weight)
        assert_array_almost_equal(
            clf_sw_lbfgs.coef_, clf_sw_n.coef_, decimal=4)
        assert_array_almost_equal(
            clf_sw_lbfgs.coef_, clf_sw_sag.coef_, decimal=4)
        assert_array_almost_equal(
            clf_sw_lbfgs.coef_, clf_sw_liblinear.coef_, decimal=4)

        # Test that passing class_weight as [1,2] is the same as
        # passing class weight = [1,1] but adjusting sample weights
        # to be 2 for all instances of class 2
        for solver in ['lbfgs', 'liblinear']:
            clf_cw_12 = LR(solver=solver, fit_intercept=False,
                           class_weight={0: 1, 1: 2}, random_state=42,
                           multi_class='ovr')
            clf_cw_12.fit(X, y)
            clf_sw_12 = LR(solver=solver, fit_intercept=False, random_state=42,
                           multi_class='ovr')
            clf_sw_12.fit(X, y, sample_weight=sample_weight)
            assert_array_almost_equal(
                clf_cw_12.coef_, clf_sw_12.coef_, decimal=4)

    # Test the above for l1 penalty and l2 penalty with dual=True.
    # since the patched liblinear code is different.
    clf_cw = LogisticRegression(
        solver="liblinear", fit_intercept=False, class_weight={0: 1, 1: 2},
        penalty="l1", tol=1e-5, random_state=42, multi_class='ovr')
    clf_cw.fit(X, y)
    clf_sw = LogisticRegression(
        solver="liblinear", fit_intercept=False, penalty="l1", tol=1e-5,
        random_state=42, multi_class='ovr')
    clf_sw.fit(X, y, sample_weight)
    assert_array_almost_equal(clf_cw.coef_, clf_sw.coef_, decimal=4)

    clf_cw = LogisticRegression(
        solver="liblinear", fit_intercept=False, class_weight={0: 1, 1: 2},
        penalty="l2", dual=True, random_state=42, multi_class='ovr')
    clf_cw.fit(X, y)
    clf_sw = LogisticRegression(
        solver="liblinear", fit_intercept=False, penalty="l2", dual=True,
        random_state=42, multi_class='ovr')
    clf_sw.fit(X, y, sample_weight)
    assert_array_almost_equal(clf_cw.coef_, clf_sw.coef_, decimal=4)


def _compute_class_weight_dictionary(y):
    # helper for returning a dictionary instead of an array
    classes = np.unique(y)
    class_weight = compute_class_weight("balanced", classes, y)
    class_weight_dict = dict(zip(classes, class_weight))
    return class_weight_dict


def test_logistic_regression_class_weights():
    # Multinomial case: remove 90% of class 0
    X = iris.data[45:, :]
    y = iris.target[45:]
    solvers = ("lbfgs", "newton-cg")
    class_weight_dict = _compute_class_weight_dictionary(y)

    for solver in solvers:
        clf1 = LogisticRegression(solver=solver, multi_class="multinomial",
                                  class_weight="balanced")
        clf2 = LogisticRegression(solver=solver, multi_class="multinomial",
                                  class_weight=class_weight_dict)
        clf1.fit(X, y)
        clf2.fit(X, y)
        assert_array_almost_equal(clf1.coef_, clf2.coef_, decimal=4)

    # Binary case: remove 90% of class 0 and 100% of class 2
    X = iris.data[45:100, :]
    y = iris.target[45:100]
    solvers = ("lbfgs", "newton-cg", "liblinear")
    class_weight_dict = _compute_class_weight_dictionary(y)

    for solver in solvers:
        clf1 = LogisticRegression(solver=solver, multi_class="ovr",
                                  class_weight="balanced")
        clf2 = LogisticRegression(solver=solver, multi_class="ovr",
                                  class_weight=class_weight_dict)
        clf1.fit(X, y)
        clf2.fit(X, y)
        assert_array_almost_equal(clf1.coef_, clf2.coef_, decimal=6)


@pytest.mark.filterwarnings('ignore: You should specify a value')  # 0.22
def test_logistic_regression_multinomial():
    # Tests for the multinomial option in logistic regression

    # Some basic attributes of Logistic Regression
    n_samples, n_features, n_classes = 50, 20, 3
    X, y = make_classification(n_samples=n_samples,
                               n_features=n_features,
                               n_informative=10,
                               n_classes=n_classes, random_state=0)

    # 'lbfgs' is used as a referenced
    solver = 'lbfgs'
    ref_i = LogisticRegression(solver=solver, multi_class='multinomial')
    ref_w = LogisticRegression(solver=solver, multi_class='multinomial',
                               fit_intercept=False)
    ref_i.fit(X, y)
    ref_w.fit(X, y)
    assert_array_equal(ref_i.coef_.shape, (n_classes, n_features))
    assert_array_equal(ref_w.coef_.shape, (n_classes, n_features))
    for solver in ['sag', 'saga', 'newton-cg']:
        clf_i = LogisticRegression(solver=solver, multi_class='multinomial',
                                   random_state=42, max_iter=2000, tol=1e-7,
                                   )
        clf_w = LogisticRegression(solver=solver, multi_class='multinomial',
                                   random_state=42, max_iter=2000, tol=1e-7,
                                   fit_intercept=False)
        clf_i.fit(X, y)
        clf_w.fit(X, y)
        assert_array_equal(clf_i.coef_.shape, (n_classes, n_features))
        assert_array_equal(clf_w.coef_.shape, (n_classes, n_features))

        # Compare solutions between lbfgs and the other solvers
        assert_almost_equal(ref_i.coef_, clf_i.coef_, decimal=3)
        assert_almost_equal(ref_w.coef_, clf_w.coef_, decimal=3)
        assert_almost_equal(ref_i.intercept_, clf_i.intercept_, decimal=3)

    # Test that the path give almost the same results. However since in this
    # case we take the average of the coefs after fitting across all the
    # folds, it need not be exactly the same.
    for solver in ['lbfgs', 'newton-cg', 'sag', 'saga']:
        clf_path = LogisticRegressionCV(solver=solver, max_iter=2000, tol=1e-6,
                                        multi_class='multinomial', Cs=[1.])
        clf_path.fit(X, y)
        assert_array_almost_equal(clf_path.coef_, ref_i.coef_, decimal=3)
        assert_almost_equal(clf_path.intercept_, ref_i.intercept_, decimal=3)


def test_multinomial_grad_hess():
    rng = np.random.RandomState(0)
    n_samples, n_features, n_classes = 100, 5, 3
    X = rng.randn(n_samples, n_features)
    w = rng.rand(n_classes, n_features)
    Y = np.zeros((n_samples, n_classes))
    ind = np.argmax(np.dot(X, w.T), axis=1)
    Y[range(0, n_samples), ind] = 1
    w = w.ravel()
    sample_weights = np.ones(X.shape[0])
    grad, hessp = _multinomial_grad_hess(w, X, Y, alpha=1.,
                                         sample_weight=sample_weights)
    # extract first column of hessian matrix
    vec = np.zeros(n_features * n_classes)
    vec[0] = 1
    hess_col = hessp(vec)

    # Estimate hessian using least squares as done in
    # test_logistic_grad_hess
    e = 1e-3
    d_x = np.linspace(-e, e, 30)
    d_grad = np.array([
        _multinomial_grad_hess(w + t * vec, X, Y, alpha=1.,
                               sample_weight=sample_weights)[0]
        for t in d_x
    ])
    d_grad -= d_grad.mean(axis=0)
    approx_hess_col = linalg.lstsq(d_x[:, np.newaxis], d_grad)[0].ravel()
    assert_array_almost_equal(hess_col, approx_hess_col)


def test_liblinear_decision_function_zero():
    # Test negative prediction when decision_function values are zero.
    # Liblinear predicts the positive class when decision_function values
    # are zero. This is a test to verify that we do not do the same.
    # See Issue: https://github.com/scikit-learn/scikit-learn/issues/3600
    # and the PR https://github.com/scikit-learn/scikit-learn/pull/3623
    X, y = make_classification(n_samples=5, n_features=5, random_state=0)
    clf = LogisticRegression(fit_intercept=False, solver='liblinear',
                             multi_class='ovr')
    clf.fit(X, y)

    # Dummy data such that the decision function becomes zero.
    X = np.zeros((5, 5))
    assert_array_equal(clf.predict(X), np.zeros(5))


@pytest.mark.filterwarnings('ignore: You should specify a value')  # 0.22
def test_liblinear_logregcv_sparse():
    # Test LogRegCV with solver='liblinear' works for sparse matrices

    X, y = make_classification(n_samples=10, n_features=5, random_state=0)
    clf = LogisticRegressionCV(solver='liblinear', multi_class='ovr')
    clf.fit(sparse.csr_matrix(X), y)


@pytest.mark.filterwarnings('ignore: Default multi_class will')  # 0.22
@pytest.mark.filterwarnings('ignore: You should specify a value')  # 0.22
def test_saga_sparse():
    # Test LogRegCV with solver='liblinear' works for sparse matrices

    X, y = make_classification(n_samples=10, n_features=5, random_state=0)
    clf = LogisticRegressionCV(solver='saga')
    clf.fit(sparse.csr_matrix(X), y)


def test_logreg_intercept_scaling():
    # Test that the right error message is thrown when intercept_scaling <= 0

    for i in [-1, 0]:
        clf = LogisticRegression(intercept_scaling=i, solver='liblinear',
                                 multi_class='ovr')
        msg = ('Intercept scaling is %r but needs to be greater than 0.'
               ' To disable fitting an intercept,'
               ' set fit_intercept=False.' % clf.intercept_scaling)
        assert_raise_message(ValueError, msg, clf.fit, X, Y1)


@pytest.mark.filterwarnings('ignore: Default solver will be changed')  # 0.22
@pytest.mark.filterwarnings('ignore: Default multi_class will')  # 0.22
def test_logreg_intercept_scaling_zero():
    # Test that intercept_scaling is ignored when fit_intercept is False

    clf = LogisticRegression(fit_intercept=False)
    clf.fit(X, Y1)
    assert_equal(clf.intercept_, 0.)


def test_logreg_l1():
    # Because liblinear penalizes the intercept and saga does not, we do not
    # fit the intercept to make it possible to compare the coefficients of
    # the two models at convergence.
    rng = np.random.RandomState(42)
    n_samples = 50
    X, y = make_classification(n_samples=n_samples, n_features=20,
                               random_state=0)
    X_noise = rng.normal(size=(n_samples, 3))
    X_constant = np.ones(shape=(n_samples, 2))
    X = np.concatenate((X, X_noise, X_constant), axis=1)
    lr_liblinear = LogisticRegression(penalty="l1", C=1.0, solver='liblinear',
                                      fit_intercept=False, multi_class='ovr',
                                      tol=1e-10)
    lr_liblinear.fit(X, y)

    lr_saga = LogisticRegression(penalty="l1", C=1.0, solver='saga',
                                 fit_intercept=False, multi_class='ovr',
                                 max_iter=1000, tol=1e-10)
    lr_saga.fit(X, y)
    assert_array_almost_equal(lr_saga.coef_, lr_liblinear.coef_)

    # Noise and constant features should be regularized to zero by the l1
    # penalty
    assert_array_almost_equal(lr_liblinear.coef_[0, -5:], np.zeros(5))
    assert_array_almost_equal(lr_saga.coef_[0, -5:], np.zeros(5))


def test_logreg_l1_sparse_data():
    # Because liblinear penalizes the intercept and saga does not, we do not
    # fit the intercept to make it possible to compare the coefficients of
    # the two models at convergence.
    rng = np.random.RandomState(42)
    n_samples = 50
    X, y = make_classification(n_samples=n_samples, n_features=20,
                               random_state=0)
    X_noise = rng.normal(scale=0.1, size=(n_samples, 3))
    X_constant = np.zeros(shape=(n_samples, 2))
    X = np.concatenate((X, X_noise, X_constant), axis=1)
    X[X < 1] = 0
    X = sparse.csr_matrix(X)

    lr_liblinear = LogisticRegression(penalty="l1", C=1.0, solver='liblinear',
                                      fit_intercept=False, multi_class='ovr',
                                      tol=1e-10)
    lr_liblinear.fit(X, y)

    lr_saga = LogisticRegression(penalty="l1", C=1.0, solver='saga',
                                 fit_intercept=False, multi_class='ovr',
                                 max_iter=1000, tol=1e-10)
    lr_saga.fit(X, y)
    assert_array_almost_equal(lr_saga.coef_, lr_liblinear.coef_)
    # Noise and constant features should be regularized to zero by the l1
    # penalty
    assert_array_almost_equal(lr_liblinear.coef_[0, -5:], np.zeros(5))
    assert_array_almost_equal(lr_saga.coef_[0, -5:], np.zeros(5))

    # Check that solving on the sparse and dense data yield the same results
    lr_saga_dense = LogisticRegression(penalty="l1", C=1.0, solver='saga',
                                       fit_intercept=False, multi_class='ovr',
                                       max_iter=1000, tol=1e-10)
    lr_saga_dense.fit(X.toarray(), y)
    assert_array_almost_equal(lr_saga.coef_, lr_saga_dense.coef_)


@pytest.mark.filterwarnings('ignore: Default multi_class will')  # 0.22
@pytest.mark.filterwarnings('ignore: You should specify a value')  # 0.22
def test_logreg_cv_penalty():
    # Test that the correct penalty is passed to the final fit.
    X, y = make_classification(n_samples=50, n_features=20, random_state=0)
    lr_cv = LogisticRegressionCV(penalty="l1", Cs=[1.0], solver='saga')
    lr_cv.fit(X, y)
    lr = LogisticRegression(penalty="l1", C=1.0, solver='saga')
    lr.fit(X, y)
    assert_equal(np.count_nonzero(lr_cv.coef_), np.count_nonzero(lr.coef_))


def test_logreg_predict_proba_multinomial():
    X, y = make_classification(n_samples=10, n_features=20, random_state=0,
                               n_classes=3, n_informative=10)

    # Predicted probabilities using the true-entropy loss should give a
    # smaller loss than those using the ovr method.
    clf_multi = LogisticRegression(multi_class="multinomial", solver="lbfgs")
    clf_multi.fit(X, y)
    clf_multi_loss = log_loss(y, clf_multi.predict_proba(X))
    clf_ovr = LogisticRegression(multi_class="ovr", solver="lbfgs")
    clf_ovr.fit(X, y)
    clf_ovr_loss = log_loss(y, clf_ovr.predict_proba(X))
    assert_greater(clf_ovr_loss, clf_multi_loss)

    # Predicted probabilities using the soft-max function should give a
    # smaller loss than those using the logistic function.
    clf_multi_loss = log_loss(y, clf_multi.predict_proba(X))
    clf_wrong_loss = log_loss(y, clf_multi._predict_proba_lr(X))
    assert_greater(clf_wrong_loss, clf_multi_loss)


def test_max_iter():
    # Test that the maximum number of iteration is reached
    X, y_bin = iris.data, iris.target.copy()
    y_bin[y_bin == 2] = 0

    solvers = ['newton-cg', 'liblinear', 'sag', 'saga', 'lbfgs']

    for max_iter in range(1, 5):
        for solver in solvers:
            for multi_class in ['ovr', 'multinomial']:
                if solver == 'liblinear' and multi_class == 'multinomial':
                    continue
                lr = LogisticRegression(max_iter=max_iter, tol=1e-15,
                                        multi_class=multi_class,
                                        random_state=0, solver=solver)
                assert_warns(ConvergenceWarning, lr.fit, X, y_bin)
                assert_equal(lr.n_iter_[0], max_iter)


@pytest.mark.parametrize('solver',
                         ['newton-cg', 'liblinear', 'sag', 'saga', 'lbfgs'])
def test_n_iter(solver):
    # Test that self.n_iter_ has the correct format.
    X, y = iris.data, iris.target
    y_bin = y.copy()
    y_bin[y_bin == 2] = 0

    n_Cs = 4
    n_cv_fold = 2

    # OvR case
    n_classes = 1 if solver == 'liblinear' else np.unique(y).shape[0]
    clf = LogisticRegression(tol=1e-2, multi_class='ovr',
                             solver=solver, C=1.,
                             random_state=42, max_iter=100)
    clf.fit(X, y)
    assert_equal(clf.n_iter_.shape, (n_classes,))

    n_classes = np.unique(y).shape[0]
    clf = LogisticRegressionCV(tol=1e-2, multi_class='ovr',
                               solver=solver, Cs=n_Cs, cv=n_cv_fold,
                               random_state=42, max_iter=100)
    clf.fit(X, y)
    assert_equal(clf.n_iter_.shape, (n_classes, n_cv_fold, n_Cs))
    clf.fit(X, y_bin)
    assert_equal(clf.n_iter_.shape, (1, n_cv_fold, n_Cs))

    # multinomial case
    n_classes = 1
    if solver in ('liblinear', 'sag', 'saga'):
        return

    clf = LogisticRegression(tol=1e-2, multi_class='multinomial',
                             solver=solver, C=1.,
                             random_state=42, max_iter=100)
    clf.fit(X, y)
    assert_equal(clf.n_iter_.shape, (n_classes,))

    clf = LogisticRegressionCV(tol=1e-2, multi_class='multinomial',
                               solver=solver, Cs=n_Cs, cv=n_cv_fold,
                               random_state=42, max_iter=100)
    clf.fit(X, y)
    assert_equal(clf.n_iter_.shape, (n_classes, n_cv_fold, n_Cs))
    clf.fit(X, y_bin)
    assert_equal(clf.n_iter_.shape, (1, n_cv_fold, n_Cs))


@pytest.mark.parametrize('solver', ('newton-cg', 'sag', 'saga', 'lbfgs'))
@pytest.mark.parametrize('warm_start', (True, False))
@pytest.mark.parametrize('fit_intercept', (True, False))
@pytest.mark.parametrize('multi_class', ['ovr', 'multinomial'])
def test_warm_start(solver, warm_start, fit_intercept, multi_class):
    # A 1-iteration second fit on same data should give almost same result
    # with warm starting, and quite different result without warm starting.
    # Warm starting does not work with liblinear solver.
    X, y = iris.data, iris.target

    clf = LogisticRegression(tol=1e-4, multi_class=multi_class,
                             warm_start=warm_start,
                             solver=solver,
                             random_state=42, max_iter=100,
                             fit_intercept=fit_intercept)
    with ignore_warnings(category=ConvergenceWarning):
        clf.fit(X, y)
        coef_1 = clf.coef_

        clf.max_iter = 1
        clf.fit(X, y)
    cum_diff = np.sum(np.abs(coef_1 - clf.coef_))
    msg = ("Warm starting issue with %s solver in %s mode "
           "with fit_intercept=%s and warm_start=%s"
           % (solver, multi_class, str(fit_intercept),
              str(warm_start)))
    if warm_start:
        assert_greater(2.0, cum_diff, msg)
    else:
        assert_greater(cum_diff, 2.0, msg)


def test_saga_vs_liblinear():
    iris = load_iris()
    X, y = iris.data, iris.target
    X = np.concatenate([X] * 10)
    y = np.concatenate([y] * 10)

    X_bin = X[y <= 1]
    y_bin = y[y <= 1] * 2 - 1

    X_sparse, y_sparse = make_classification(n_samples=50, n_features=20,
                                             random_state=0)
    X_sparse = sparse.csr_matrix(X_sparse)

    for (X, y) in ((X_bin, y_bin), (X_sparse, y_sparse)):
        for penalty in ['l1', 'l2']:
            n_samples = X.shape[0]
            # alpha=1e-3 is time consuming
            for alpha in np.logspace(-1, 1, 3):
                saga = LogisticRegression(
                    C=1. / (n_samples * alpha),
                    solver='saga',
                    multi_class='ovr',
                    max_iter=200,
                    fit_intercept=False,
                    penalty=penalty, random_state=0, tol=1e-24)

                liblinear = LogisticRegression(
                    C=1. / (n_samples * alpha),
                    solver='liblinear',
                    multi_class='ovr',
                    max_iter=200,
                    fit_intercept=False,
                    penalty=penalty, random_state=0, tol=1e-24)

                saga.fit(X, y)
                liblinear.fit(X, y)
                # Convergence for alpha=1e-3 is very slow
                assert_array_almost_equal(saga.coef_, liblinear.coef_, 3)


@pytest.mark.parametrize('multi_class', ['ovr', 'multinomial'])
def test_dtype_match(multi_class):
    # Test that np.float32 input data is not cast to np.float64 when possible

    X_32 = np.array(X).astype(np.float32)
    y_32 = np.array(Y1).astype(np.float32)
    X_64 = np.array(X).astype(np.float64)
    y_64 = np.array(Y1).astype(np.float64)
    X_sparse_32 = sp.csr_matrix(X, dtype=np.float32)

    solver = 'newton-cg'

    # Check type consistency
    lr_32 = LogisticRegression(solver=solver, multi_class=multi_class,
                               random_state=42)
    lr_32.fit(X_32, y_32)
    assert_equal(lr_32.coef_.dtype, X_32.dtype)

    # check consistency with sparsity
    lr_32_sparse = LogisticRegression(solver=solver,
                                      multi_class=multi_class,
                                      random_state=42)
    lr_32_sparse.fit(X_sparse_32, y_32)
    assert_equal(lr_32_sparse.coef_.dtype, X_sparse_32.dtype)

    # Check accuracy consistency
    lr_64 = LogisticRegression(solver=solver, multi_class=multi_class,
                               random_state=42)
    lr_64.fit(X_64, y_64)
    assert_equal(lr_64.coef_.dtype, X_64.dtype)

    rtol = 1e-6
    if os.name == 'nt' and _IS_32BIT:
        # FIXME
        rtol = 1e-2

    assert_allclose(lr_32.coef_, lr_64.coef_.astype(np.float32), rtol=rtol)


def test_warm_start_converge_LR():
    # Test to see that the logistic regression converges on warm start,
    # with multi_class='multinomial'. Non-regressive test for #10836

    rng = np.random.RandomState(0)
    X = np.concatenate((rng.randn(100, 2) + [1, 1], rng.randn(100, 2)))
    y = np.array([1] * 100 + [-1] * 100)
    lr_no_ws = LogisticRegression(multi_class='multinomial',
                                  solver='sag', warm_start=False,
                                  random_state=0)
    lr_ws = LogisticRegression(multi_class='multinomial',
                               solver='sag', warm_start=True,
                               random_state=0)

    lr_no_ws_loss = log_loss(y, lr_no_ws.fit(X, y).predict_proba(X))
    for i in range(5):
        lr_ws.fit(X, y)
    lr_ws_loss = log_loss(y, lr_ws.predict_proba(X))
    assert_allclose(lr_no_ws_loss, lr_ws_loss, rtol=1e-5)


def test_logistic_regression_path_coefs_multinomial():
    # Make sure that the returned coefs by logistic_regression_path when
    # multi_class='multinomial' don't override each other (used to be a
    # bug).
    X, y = make_classification(n_samples=200, n_classes=3, n_informative=2,
                               n_redundant=0, n_clusters_per_class=1,
                               random_state=0, n_features=2)
    Cs = [.00001, 1, 10000]
    coefs, _, _ = logistic_regression_path(X, y, penalty='l1', Cs=Cs,
                                           solver='saga', random_state=0,
                                           multi_class='multinomial')

    with pytest.raises(AssertionError):
        assert_array_almost_equal(coefs[0], coefs[1], decimal=1)
    with pytest.raises(AssertionError):
        assert_array_almost_equal(coefs[0], coefs[2], decimal=1)
    with pytest.raises(AssertionError):
        assert_array_almost_equal(coefs[1], coefs[2], decimal=1)


@pytest.mark.parametrize('est', [LogisticRegression(random_state=0),
                                 LogisticRegressionCV(random_state=0, cv=3),
                                 ])
@pytest.mark.parametrize('solver', ['liblinear', 'lbfgs', 'newton-cg', 'sag',
                                    'saga'])
def test_logistic_regression_multi_class_auto(est, solver):
    # check multi_class='auto' => multi_class='ovr' iff binary y or liblinear

    def fit(X, y, **kw):
        return clone(est).set_params(**kw).fit(X, y)

    X = iris.data[::10]
    X2 = iris.data[1::10]
    y_multi = iris.target[::10]
    y_bin = y_multi == 0
    est_auto_bin = fit(X, y_bin, multi_class='auto', solver=solver)
    est_ovr_bin = fit(X, y_bin, multi_class='ovr', solver=solver)
    assert np.allclose(est_auto_bin.coef_, est_ovr_bin.coef_)
    assert np.allclose(est_auto_bin.predict_proba(X2),
                       est_ovr_bin.predict_proba(X2))

    est_auto_multi = fit(X, y_multi, multi_class='auto', solver=solver)
    if solver == 'liblinear':
        est_ovr_multi = fit(X, y_multi, multi_class='ovr', solver=solver)
        assert np.allclose(est_auto_multi.coef_, est_ovr_multi.coef_)
        assert np.allclose(est_auto_multi.predict_proba(X2),
                           est_ovr_multi.predict_proba(X2))
    else:
        est_multi_multi = fit(X, y_multi, multi_class='multinomial',
                              solver=solver)
        if sys.platform == 'darwin' and solver == 'lbfgs':
            pytest.xfail('Issue #11924: LogisticRegressionCV(solver="lbfgs", '
                         'multi_class="multinomial") is nondterministic on '
                         'MacOS.')  # pragma: no cover
        assert np.allclose(est_auto_multi.coef_, est_multi_multi.coef_)
        assert np.allclose(est_auto_multi.predict_proba(X2),
                           est_multi_multi.predict_proba(X2))

        # Make sure multi_class='ovr' is distinct from ='multinomial'
        assert not np.allclose(est_auto_bin.coef_,
                               fit(X, y_bin, multi_class='multinomial',
                                   solver=solver).coef_)
        assert not np.allclose(est_auto_bin.coef_,
                               fit(X, y_multi, multi_class='multinomial',
                                   solver=solver).coef_)
