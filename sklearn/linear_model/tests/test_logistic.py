import itertools
import os
import warnings
from functools import partial

import numpy as np
import pytest
from numpy.testing import (
    assert_allclose,
    assert_almost_equal,
    assert_array_almost_equal,
    assert_array_equal,
)
from scipy import sparse
from scipy.linalg import LinAlgWarning, svd

from sklearn import config_context
from sklearn._loss import HalfMultinomialLoss
from sklearn.base import clone
from sklearn.datasets import load_iris, make_classification, make_low_rank_matrix
from sklearn.exceptions import ConvergenceWarning
from sklearn.linear_model import SGDClassifier
from sklearn.linear_model._logistic import (
    LogisticRegression as LogisticRegressionDefault,
)
from sklearn.linear_model._logistic import (
    LogisticRegressionCV as LogisticRegressionCVDefault,
)
from sklearn.linear_model._logistic import (
    _log_reg_scoring_path,
    _logistic_regression_path,
)
from sklearn.metrics import get_scorer, log_loss
from sklearn.model_selection import (
    GridSearchCV,
    LeaveOneGroupOut,
    StratifiedKFold,
    cross_val_score,
    train_test_split,
)
from sklearn.multiclass import OneVsRestClassifier
from sklearn.preprocessing import LabelEncoder, StandardScaler, scale
from sklearn.svm import l1_min_c
from sklearn.utils import compute_class_weight, shuffle
from sklearn.utils._testing import ignore_warnings, skip_if_no_parallel
from sklearn.utils.fixes import _IS_32BIT, COO_CONTAINERS, CSR_CONTAINERS

pytestmark = pytest.mark.filterwarnings(
    "error::sklearn.exceptions.ConvergenceWarning:sklearn.*"
)
# Fixing random_state helps prevent ConvergenceWarnings
LogisticRegression = partial(LogisticRegressionDefault, random_state=0)
LogisticRegressionCV = partial(LogisticRegressionCVDefault, random_state=0)


SOLVERS = ("lbfgs", "liblinear", "newton-cg", "newton-cholesky", "sag", "saga")
X = [[-1, 0], [0, 1], [1, 1]]
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

    assert predicted.shape == (n_samples,)
    assert_array_equal(predicted, y)

    probabilities = clf.predict_proba(X)
    assert probabilities.shape == (n_samples, n_classes)
    assert_array_almost_equal(probabilities.sum(axis=1), np.ones(n_samples))
    assert_array_equal(probabilities.argmax(axis=1), y)


@pytest.mark.parametrize("csr_container", CSR_CONTAINERS)
def test_predict_2_classes(csr_container):
    # Simple sanity check on a 2 classes dataset
    # Make sure it predicts the correct result on simple datasets.
    check_predictions(LogisticRegression(random_state=0), X, Y1)
    check_predictions(LogisticRegression(random_state=0), csr_container(X), Y1)

    check_predictions(LogisticRegression(C=100, random_state=0), X, Y1)
    check_predictions(LogisticRegression(C=100, random_state=0), csr_container(X), Y1)

    check_predictions(LogisticRegression(fit_intercept=False, random_state=0), X, Y1)
    check_predictions(
        LogisticRegression(fit_intercept=False, random_state=0), csr_container(X), Y1
    )


def test_logistic_cv_mock_scorer():
    class MockScorer:
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
    X, y = make_classification(random_state=0)
    lr.fit(X, y)

    # Cs[2] has the highest score (0.8) from MockScorer
    assert lr.C_[0] == Cs[2]

    # scorer called 8 times (cv*len(Cs))
    assert mock_scorer.calls == cv * len(Cs)

    # reset mock_scorer
    mock_scorer.calls = 0
    custom_score = lr.score(X, lr.predict(X))

    assert custom_score == mock_scorer.scores[0]
    assert mock_scorer.calls == 1


@skip_if_no_parallel
def test_lr_liblinear_warning():
    X, y = make_classification(random_state=0)

    lr = LogisticRegression(solver="liblinear", n_jobs=2)
    warning_message = (
        "'n_jobs' > 1 does not have any effect when"
        " 'solver' is set to 'liblinear'. Got 'n_jobs'"
        " = 2."
    )
    with pytest.warns(UserWarning, match=warning_message):
        lr.fit(X, y)


@pytest.mark.parametrize("csr_container", CSR_CONTAINERS)
def test_predict_3_classes(csr_container):
    check_predictions(LogisticRegression(C=10), X, Y2)
    check_predictions(LogisticRegression(C=10), csr_container(X), Y2)


# TODO(1.8): remove filterwarnings after the deprecation of multi_class
@pytest.mark.filterwarnings("ignore:.*'multi_class' was deprecated.*:FutureWarning")
@pytest.mark.filterwarnings(
    "ignore:.*'liblinear' solver for multiclass classification is deprecated.*"
)
@pytest.mark.parametrize(
    "clf",
    [
        LogisticRegression(C=len(iris.data), solver="liblinear", multi_class="ovr"),
        LogisticRegression(C=len(iris.data), solver="lbfgs"),
        LogisticRegression(C=len(iris.data), solver="newton-cg"),
        LogisticRegression(
            C=len(iris.data), solver="sag", tol=1e-2, multi_class="ovr", random_state=42
        ),
        LogisticRegression(
            C=len(iris.data),
            solver="saga",
            tol=1e-2,
            multi_class="ovr",
            random_state=42,
        ),
        LogisticRegression(C=len(iris.data), solver="newton-cholesky"),
    ],
)
def test_predict_iris(clf):
    """Test logistic regression with the iris dataset.

    Test that both multinomial and OvR solvers handle multiclass data correctly and
    give good accuracy score (>0.95) for the training data.
    """
    n_samples, n_features = iris.data.shape
    target = iris.target_names[iris.target]

    if clf.solver == "lbfgs":
        # lbfgs has convergence issues on the iris data with its default max_iter=100
        with warnings.catch_warnings():
            warnings.simplefilter("ignore", ConvergenceWarning)
            clf.fit(iris.data, target)
    else:
        clf.fit(iris.data, target)
    assert_array_equal(np.unique(target), clf.classes_)

    pred = clf.predict(iris.data)
    assert np.mean(pred == target) > 0.95

    probabilities = clf.predict_proba(iris.data)
    assert_allclose(probabilities.sum(axis=1), np.ones(n_samples))

    pred = iris.target_names[probabilities.argmax(axis=1)]
    assert np.mean(pred == target) > 0.95


# TODO(1.8): remove filterwarnings after the deprecation of multi_class
@pytest.mark.filterwarnings("ignore:.*'multi_class' was deprecated.*:FutureWarning")
@pytest.mark.parametrize("LR", [LogisticRegression, LogisticRegressionCV])
def test_check_solver_option(LR):
    X, y = iris.data, iris.target

    # only 'liblinear' solver
    for solver in ["liblinear"]:
        msg = f"Solver {solver} does not support a multinomial backend."
        lr = LR(solver=solver, multi_class="multinomial")
        with pytest.raises(ValueError, match=msg):
            lr.fit(X, y)

    # all solvers except 'liblinear' and 'saga'
    for solver in ["lbfgs", "newton-cg", "newton-cholesky", "sag"]:
        msg = "Solver %s supports only 'l2' or None penalties," % solver
        lr = LR(solver=solver, penalty="l1", multi_class="ovr")
        with pytest.raises(ValueError, match=msg):
            lr.fit(X, y)
    for solver in ["lbfgs", "newton-cg", "newton-cholesky", "sag", "saga"]:
        msg = "Solver %s supports only dual=False, got dual=True" % solver
        lr = LR(solver=solver, dual=True, multi_class="ovr")
        with pytest.raises(ValueError, match=msg):
            lr.fit(X, y)

    # only saga supports elasticnet. We only test for liblinear because the
    # error is raised before for the other solvers (solver %s supports only l2
    # penalties)
    for solver in ["liblinear"]:
        msg = f"Only 'saga' solver supports elasticnet penalty, got solver={solver}."
        lr = LR(solver=solver, penalty="elasticnet")
        with pytest.raises(ValueError, match=msg):
            lr.fit(X, y)

    # liblinear does not support penalty='none'
    # (LogisticRegressionCV does not supports penalty='none' at all)
    if LR is LogisticRegression:
        msg = "penalty=None is not supported for the liblinear solver"
        lr = LR(penalty=None, solver="liblinear")
        with pytest.raises(ValueError, match=msg):
            lr.fit(X, y)


@pytest.mark.parametrize("LR", [LogisticRegression, LogisticRegressionCV])
def test_elasticnet_l1_ratio_err_helpful(LR):
    # Check that an informative error message is raised when penalty="elasticnet"
    # but l1_ratio is not specified.
    model = LR(penalty="elasticnet", solver="saga")
    with pytest.raises(ValueError, match=r".*l1_ratio.*"):
        model.fit(np.array([[1, 2], [3, 4]]), np.array([0, 1]))


# TODO(1.8): remove whole test with deprecation of multi_class
@pytest.mark.filterwarnings("ignore:.*'multi_class' was deprecated.*:FutureWarning")
@pytest.mark.parametrize("solver", ["lbfgs", "newton-cg", "sag", "saga"])
def test_multinomial_binary(solver):
    # Test multinomial LR on a binary problem.
    target = (iris.target > 0).astype(np.intp)
    target = np.array(["setosa", "not-setosa"])[target]

    clf = LogisticRegression(
        solver=solver, multi_class="multinomial", random_state=42, max_iter=2000
    )
    clf.fit(iris.data, target)

    assert clf.coef_.shape == (1, iris.data.shape[1])
    assert clf.intercept_.shape == (1,)
    assert_array_equal(clf.predict(iris.data), target)

    mlr = LogisticRegression(
        solver=solver, multi_class="multinomial", random_state=42, fit_intercept=False
    )
    mlr.fit(iris.data, target)
    pred = clf.classes_[np.argmax(clf.predict_log_proba(iris.data), axis=1)]
    assert np.mean(pred == target) > 0.9


# TODO(1.8): remove filterwarnings after the deprecation of multi_class
# Maybe even remove this whole test as correctness of multinomial loss is tested
# elsewhere.
@pytest.mark.filterwarnings("ignore:.*'multi_class' was deprecated.*:FutureWarning")
def test_multinomial_binary_probabilities(global_random_seed):
    # Test multinomial LR gives expected probabilities based on the
    # decision function, for a binary problem.
    X, y = make_classification(random_state=global_random_seed)
    clf = LogisticRegression(
        multi_class="multinomial",
        solver="saga",
        tol=1e-3,
        random_state=global_random_seed,
    )
    clf.fit(X, y)

    decision = clf.decision_function(X)
    proba = clf.predict_proba(X)

    expected_proba_class_1 = np.exp(decision) / (np.exp(decision) + np.exp(-decision))
    expected_proba = np.c_[1 - expected_proba_class_1, expected_proba_class_1]

    assert_almost_equal(proba, expected_proba)


@pytest.mark.parametrize("coo_container", COO_CONTAINERS)
def test_sparsify(coo_container):
    # Test sparsify and densify members.
    n_samples, n_features = iris.data.shape
    target = iris.target_names[iris.target]
    X = scale(iris.data)
    clf = LogisticRegression(random_state=0).fit(X, target)

    pred_d_d = clf.decision_function(X)

    clf.sparsify()
    assert sparse.issparse(clf.coef_)
    pred_s_d = clf.decision_function(X)

    sp_data = coo_container(X)
    pred_s_s = clf.decision_function(sp_data)

    clf.densify()
    pred_d_s = clf.decision_function(sp_data)

    assert_array_almost_equal(pred_d_d, pred_s_d)
    assert_array_almost_equal(pred_d_d, pred_s_s)
    assert_array_almost_equal(pred_d_d, pred_d_s)


def test_inconsistent_input():
    # Test that an exception is raised on inconsistent input
    rng = np.random.RandomState(0)
    X_ = rng.random_sample((5, 10))
    y_ = np.ones(X_.shape[0])
    y_[0] = 0

    clf = LogisticRegression(random_state=0)

    # Wrong dimensions for training data
    y_wrong = y_[:-1]

    with pytest.raises(ValueError):
        clf.fit(X, y_wrong)

    # Wrong dimensions for test data
    with pytest.raises(ValueError):
        clf.fit(X_, y_).predict(rng.random_sample((3, 12)))


def test_write_parameters():
    # Test that we can write to coef_ and intercept_
    clf = LogisticRegression(random_state=0)
    clf.fit(X, Y1)
    clf.coef_[:] = 0
    clf.intercept_[:] = 0
    assert_array_almost_equal(clf.decision_function(X), 0)


def test_nan():
    # Test proper NaN handling.
    # Regression test for Issue #252: fit used to go into an infinite loop.
    Xnan = np.array(X, dtype=np.float64)
    Xnan[0, 1] = np.nan
    logistic = LogisticRegression(random_state=0)

    with pytest.raises(ValueError):
        logistic.fit(Xnan, Y1)


def test_consistency_path():
    # Test that the path algorithm is consistent
    rng = np.random.RandomState(0)
    X = np.concatenate((rng.randn(100, 2) + [1, 1], rng.randn(100, 2)))
    y = [1] * 100 + [-1] * 100
    Cs = np.logspace(0, 4, 10)

    f = ignore_warnings
    # can't test with fit_intercept=True since LIBLINEAR
    # penalizes the intercept
    for solver in ["sag", "saga"]:
        coefs, Cs, _ = f(_logistic_regression_path)(
            X,
            y,
            Cs=Cs,
            fit_intercept=False,
            tol=1e-5,
            solver=solver,
            max_iter=1000,
            random_state=0,
        )
        for i, C in enumerate(Cs):
            lr = LogisticRegression(
                C=C,
                fit_intercept=False,
                tol=1e-5,
                solver=solver,
                random_state=0,
                max_iter=1000,
            )
            lr.fit(X, y)
            lr_coef = lr.coef_.ravel()
            assert_array_almost_equal(
                lr_coef, coefs[i], decimal=4, err_msg="with solver = %s" % solver
            )

    # test for fit_intercept=True
    for solver in ("lbfgs", "newton-cg", "newton-cholesky", "liblinear", "sag", "saga"):
        Cs = [1e3]
        coefs, Cs, _ = f(_logistic_regression_path)(
            X,
            y,
            Cs=Cs,
            tol=1e-6,
            solver=solver,
            intercept_scaling=10000.0,
            random_state=0,
        )
        lr = LogisticRegression(
            C=Cs[0],
            tol=1e-6,
            intercept_scaling=10000.0,
            random_state=0,
            solver=solver,
        )
        lr.fit(X, y)
        lr_coef = np.concatenate([lr.coef_.ravel(), lr.intercept_])
        assert_array_almost_equal(
            lr_coef, coefs[0], decimal=4, err_msg="with solver = %s" % solver
        )


def test_logistic_regression_path_convergence_fail():
    rng = np.random.RandomState(0)
    X = np.concatenate((rng.randn(100, 2) + [1, 1], rng.randn(100, 2)))
    y = [1] * 100 + [-1] * 100
    Cs = [1e3]

    # Check that the convergence message points to both a model agnostic
    # advice (scaling the data) and to the logistic regression specific
    # documentation that includes hints on the solver configuration.
    with pytest.warns(ConvergenceWarning) as record:
        _logistic_regression_path(
            X, y, Cs=Cs, tol=0.0, max_iter=1, random_state=0, verbose=0
        )

    assert len(record) == 1
    warn_msg = record[0].message.args[0]
    assert "lbfgs failed to converge after 1 iteration(s)" in warn_msg
    assert "Increase the number of iterations" in warn_msg
    assert "scale the data" in warn_msg
    assert "linear_model.html#logistic-regression" in warn_msg


def test_liblinear_dual_random_state():
    # random_state is relevant for liblinear solver only if dual=True
    X, y = make_classification(n_samples=20, random_state=0)
    lr1 = LogisticRegression(
        random_state=0,
        dual=True,
        tol=1e-3,
        solver="liblinear",
    )
    lr1.fit(X, y)
    lr2 = LogisticRegression(
        random_state=0,
        dual=True,
        tol=1e-3,
        solver="liblinear",
    )
    lr2.fit(X, y)
    lr3 = LogisticRegression(
        random_state=8,
        dual=True,
        tol=1e-3,
        solver="liblinear",
    )
    lr3.fit(X, y)

    # same result for same random state
    assert_array_almost_equal(lr1.coef_, lr2.coef_)
    # different results for different random states
    msg = "Arrays are not almost equal to 6 decimals"
    with pytest.raises(AssertionError, match=msg):
        assert_array_almost_equal(lr1.coef_, lr3.coef_)


def test_logistic_cv():
    # test for LogisticRegressionCV object
    n_samples, n_features = 50, 5
    rng = np.random.RandomState(0)
    X_ref = rng.randn(n_samples, n_features)
    y = np.sign(X_ref.dot(5 * rng.randn(n_features)))
    X_ref -= X_ref.mean()
    X_ref /= X_ref.std()
    lr_cv = LogisticRegressionCV(
        Cs=[1.0], fit_intercept=False, solver="liblinear", cv=3
    )
    lr_cv.fit(X_ref, y)
    lr = LogisticRegression(C=1.0, fit_intercept=False, solver="liblinear")
    lr.fit(X_ref, y)
    assert_array_almost_equal(lr.coef_, lr_cv.coef_)

    assert_array_equal(lr_cv.coef_.shape, (1, n_features))
    assert_array_equal(lr_cv.classes_, [-1, 1])
    assert len(lr_cv.classes_) == 2

    coefs_paths = np.asarray(list(lr_cv.coefs_paths_.values()))
    assert_array_equal(coefs_paths.shape, (1, 3, 1, n_features))
    assert_array_equal(lr_cv.Cs_.shape, (1,))
    scores = np.asarray(list(lr_cv.scores_.values()))
    assert_array_equal(scores.shape, (1, 3, 1))


@pytest.mark.parametrize(
    "scoring, multiclass_agg_list",
    [
        ("accuracy", [""]),
        ("precision", ["_macro", "_weighted"]),
        # no need to test for micro averaging because it
        # is the same as accuracy for f1, precision,
        # and recall (see https://github.com/
        # scikit-learn/scikit-learn/pull/
        # 11578#discussion_r203250062)
        ("f1", ["_macro", "_weighted"]),
        ("neg_log_loss", [""]),
        ("recall", ["_macro", "_weighted"]),
    ],
)
def test_logistic_cv_multinomial_score(scoring, multiclass_agg_list):
    # test that LogisticRegressionCV uses the right score to compute its
    # cross-validation scores when using a multinomial scoring
    # see https://github.com/scikit-learn/scikit-learn/issues/8720
    X, y = make_classification(
        n_samples=100, random_state=0, n_classes=3, n_informative=6
    )
    train, test = np.arange(80), np.arange(80, 100)
    lr = LogisticRegression(C=1.0)
    # we use lbfgs to support multinomial
    params = lr.get_params()
    # we store the params to set them further in _log_reg_scoring_path
    for key in ["C", "n_jobs", "warm_start"]:
        del params[key]
    lr.fit(X[train], y[train])
    for averaging in multiclass_agg_list:
        scorer = get_scorer(scoring + averaging)
        assert_array_almost_equal(
            _log_reg_scoring_path(
                X,
                y,
                train,
                test,
                Cs=[1.0],
                scoring=scorer,
                pos_class=None,
                max_squared_sum=None,
                sample_weight=None,
                score_params=None,
                **(params | {"multi_class": "multinomial"}),
            )[2][0],
            scorer(lr, X[test], y[test]),
        )


def test_multinomial_logistic_regression_string_inputs():
    # Test with string labels for LogisticRegression(CV)
    n_samples, n_features, n_classes = 50, 5, 3
    X_ref, y = make_classification(
        n_samples=n_samples,
        n_features=n_features,
        n_classes=n_classes,
        n_informative=3,
        random_state=0,
    )
    y_str = LabelEncoder().fit(["bar", "baz", "foo"]).inverse_transform(y)
    # For numerical labels, let y values be taken from set (-1, 0, 1)
    y = np.array(y) - 1
    # Test for string labels
    lr = LogisticRegression()
    lr_cv = LogisticRegressionCV(Cs=3)
    lr_str = LogisticRegression()
    lr_cv_str = LogisticRegressionCV(Cs=3)

    lr.fit(X_ref, y)
    lr_cv.fit(X_ref, y)
    lr_str.fit(X_ref, y_str)
    lr_cv_str.fit(X_ref, y_str)

    assert_array_almost_equal(lr.coef_, lr_str.coef_)
    assert sorted(lr_str.classes_) == ["bar", "baz", "foo"]
    assert_array_almost_equal(lr_cv.coef_, lr_cv_str.coef_)
    assert sorted(lr_str.classes_) == ["bar", "baz", "foo"]
    assert sorted(lr_cv_str.classes_) == ["bar", "baz", "foo"]

    # The predictions should be in original labels
    assert sorted(np.unique(lr_str.predict(X_ref))) == ["bar", "baz", "foo"]
    assert sorted(np.unique(lr_cv_str.predict(X_ref))) == ["bar", "baz", "foo"]

    # Make sure class weights can be given with string labels
    lr_cv_str = LogisticRegression(class_weight={"bar": 1, "baz": 2, "foo": 0}).fit(
        X_ref, y_str
    )
    assert sorted(np.unique(lr_cv_str.predict(X_ref))) == ["bar", "baz"]


@pytest.mark.parametrize("csr_container", CSR_CONTAINERS)
def test_logistic_cv_sparse(csr_container):
    X, y = make_classification(n_samples=50, n_features=5, random_state=0)
    X[X < 1.0] = 0.0
    csr = csr_container(X)

    clf = LogisticRegressionCV()
    clf.fit(X, y)
    clfs = LogisticRegressionCV()
    clfs.fit(csr, y)
    assert_array_almost_equal(clfs.coef_, clf.coef_)
    assert_array_almost_equal(clfs.intercept_, clf.intercept_)
    assert clfs.C_ == clf.C_


# TODO(1.8): remove filterwarnings after the deprecation of multi_class
# Best remove this whole test.
@pytest.mark.filterwarnings("ignore:.*'multi_class' was deprecated.*:FutureWarning")
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
    clf = LogisticRegressionCV(cv=precomputed_folds, multi_class="ovr")
    clf.fit(train, target)

    # Conflate classes 0 and 1 and train clf1 on this modified dataset
    clf1 = LogisticRegressionCV(cv=precomputed_folds, multi_class="ovr")
    target_copy = target.copy()
    target_copy[target_copy == 0] = 1
    clf1.fit(train, target_copy)

    # Ensure that what OvR learns for class2 is same regardless of whether
    # classes 0 and 1 are separated or not
    assert_allclose(clf.scores_[2], clf1.scores_[2])
    assert_allclose(clf.intercept_[2:], clf1.intercept_)
    assert_allclose(clf.coef_[2][np.newaxis, :], clf1.coef_)

    # Test the shape of various attributes.
    assert clf.coef_.shape == (3, n_features)
    assert_array_equal(clf.classes_, [0, 1, 2])
    coefs_paths = np.asarray(list(clf.coefs_paths_.values()))
    assert coefs_paths.shape == (3, n_cv, 10, n_features + 1)
    assert clf.Cs_.shape == (10,)
    scores = np.asarray(list(clf.scores_.values()))
    assert scores.shape == (3, n_cv, 10)

    # Test that for the iris data multinomial gives a better accuracy than OvR
    for solver in ["lbfgs", "newton-cg", "sag", "saga"]:
        max_iter = 500 if solver in ["sag", "saga"] else 30
        clf_multi = LogisticRegressionCV(
            solver=solver,
            max_iter=max_iter,
            random_state=42,
            tol=1e-3 if solver in ["sag", "saga"] else 1e-2,
            cv=2,
        )
        if solver == "lbfgs":
            # lbfgs requires scaling to avoid convergence warnings
            train = scale(train)

        clf_multi.fit(train, target)
        multi_score = clf_multi.score(train, target)
        ovr_score = clf.score(train, target)
        assert multi_score > ovr_score

        # Test attributes of LogisticRegressionCV
        assert clf.coef_.shape == clf_multi.coef_.shape
        assert_array_equal(clf_multi.classes_, [0, 1, 2])
        coefs_paths = np.asarray(list(clf_multi.coefs_paths_.values()))
        assert coefs_paths.shape == (3, n_cv, 10, n_features + 1)
        assert clf_multi.Cs_.shape == (10,)
        scores = np.asarray(list(clf_multi.scores_.values()))
        assert scores.shape == (3, n_cv, 10)


def test_logistic_regression_solvers():
    """Test solvers converge to the same result."""
    X, y = make_classification(n_features=10, n_informative=5, random_state=0)

    params = dict(fit_intercept=False, random_state=42)

    regressors = {
        solver: LogisticRegression(solver=solver, **params).fit(X, y)
        for solver in SOLVERS
    }

    for solver_1, solver_2 in itertools.combinations(regressors, r=2):
        assert_array_almost_equal(
            regressors[solver_1].coef_, regressors[solver_2].coef_, decimal=3
        )


# TODO(1.8): remove filterwarnings after the deprecation of multi_class
@pytest.mark.filterwarnings("ignore:.*'multi_class' was deprecated.*:FutureWarning")
@pytest.mark.parametrize("fit_intercept", [False, True])
def test_logistic_regression_solvers_multiclass(fit_intercept):
    """Test solvers converge to the same result for multiclass problems."""
    X, y = make_classification(
        n_samples=20, n_features=20, n_informative=10, n_classes=3, random_state=0
    )
    tol = 1e-8
    params = dict(fit_intercept=fit_intercept, tol=tol, random_state=42)

    # Override max iteration count for specific solvers to allow for
    # proper convergence.
    solver_max_iter = {"lbfgs": 200, "sag": 10_000, "saga": 10_000}

    regressors = {
        solver: LogisticRegression(
            solver=solver, max_iter=solver_max_iter.get(solver, 100), **params
        ).fit(X, y)
        for solver in set(SOLVERS) - set(["liblinear"])
    }

    for solver_1, solver_2 in itertools.combinations(regressors, r=2):
        assert_allclose(
            regressors[solver_1].coef_,
            regressors[solver_2].coef_,
            rtol=5e-3 if (solver_1 == "saga" or solver_2 == "saga") else 1e-3,
            err_msg=f"{solver_1} vs {solver_2}",
        )
        if fit_intercept:
            assert_allclose(
                regressors[solver_1].intercept_,
                regressors[solver_2].intercept_,
                rtol=5e-3 if (solver_1 == "saga" or solver_2 == "saga") else 1e-3,
                err_msg=f"{solver_1} vs {solver_2}",
            )


@pytest.mark.parametrize("fit_intercept", [False, True])
def test_logistic_regression_solvers_multiclass_unpenalized(
    fit_intercept, global_random_seed
):
    """Test and compare solver results for unpenalized multinomial multiclass."""
    # We want to avoid perfect separation.
    n_samples, n_features, n_classes = 100, 4, 3
    rng = np.random.RandomState(global_random_seed)
    X = make_low_rank_matrix(
        n_samples=n_samples,
        n_features=n_features + fit_intercept,
        effective_rank=n_features + fit_intercept,
        tail_strength=0.1,
        random_state=rng,
    )
    if fit_intercept:
        X[:, -1] = 1
    U, s, Vt = svd(X)
    assert np.all(s > 1e-3)  # to be sure that X is not singular
    assert np.max(s) / np.min(s) < 100  # condition number of X
    if fit_intercept:
        X = X[:, :-1]
    coef = rng.uniform(low=1, high=3, size=n_features * n_classes)
    coef = coef.reshape(n_classes, n_features)
    intercept = rng.uniform(low=-1, high=1, size=n_classes) * fit_intercept
    raw_prediction = X @ coef.T + intercept

    loss = HalfMultinomialLoss(n_classes=n_classes)
    proba = loss.link.inverse(raw_prediction)
    # Only newer numpy version (1.22) support more dimensions on pvals.
    y = np.zeros(n_samples)
    for i in range(n_samples):
        y[i] = np.argwhere(rng.multinomial(n=1, pvals=proba[i, :]))[0, 0]

    tol = 1e-9
    params = dict(fit_intercept=fit_intercept, random_state=42)
    solver_max_iter = {"lbfgs": 200, "sag": 10_000, "saga": 10_000}
    solver_tol = {"sag": 1e-8, "saga": 1e-8}
    regressors = {
        solver: LogisticRegression(
            C=np.inf,
            solver=solver,
            tol=solver_tol.get(solver, tol),
            max_iter=solver_max_iter.get(solver, 100),
            **params,
        ).fit(X, y)
        for solver in set(SOLVERS) - set(["liblinear"])
    }
    for solver in regressors.keys():
        # See the docstring of test_multinomial_identifiability_on_iris for reference.
        assert_allclose(
            regressors[solver].coef_.sum(axis=0), 0, atol=1e-10, err_msg=solver
        )

    for solver_1, solver_2 in itertools.combinations(regressors, r=2):
        assert_allclose(
            regressors[solver_1].coef_,
            regressors[solver_2].coef_,
            rtol=5e-3 if (solver_1 == "saga" or solver_2 == "saga") else 2e-3,
            err_msg=f"{solver_1} vs {solver_2}",
        )
        if fit_intercept:
            assert_allclose(
                regressors[solver_1].intercept_,
                regressors[solver_2].intercept_,
                rtol=5e-3 if (solver_1 == "saga" or solver_2 == "saga") else 1e-3,
                err_msg=f"{solver_1} vs {solver_2}",
            )


@pytest.mark.parametrize("weight", [{0: 0.1, 1: 0.2}, {0: 0.1, 1: 0.2, 2: 0.5}])
@pytest.mark.parametrize("class_weight", ["weight", "balanced"])
def test_logistic_regressioncv_class_weights(weight, class_weight, global_random_seed):
    """Test class_weight for LogisticRegressionCV."""
    n_classes = len(weight)
    if class_weight == "weight":
        class_weight = weight

    X, y = make_classification(
        n_samples=30,
        n_features=3,
        n_repeated=0,
        n_informative=3,
        n_redundant=0,
        n_classes=n_classes,
        random_state=global_random_seed,
    )
    params = dict(
        Cs=1,
        fit_intercept=False,
        class_weight=class_weight,
        tol=1e-8,
    )
    clf_lbfgs = LogisticRegressionCV(solver="lbfgs", **params)

    # XXX: lbfgs' line search can fail and cause a ConvergenceWarning for some
    # 10% of the random seeds, but only on specific platforms (in particular
    # when using Atlas BLAS/LAPACK implementation). Doubling the maxls internal
    # parameter of the solver does not help. However this lack of proper
    # convergence does not seem to prevent the assertion to pass, so we ignore
    # the warning for now.
    # See: https://github.com/scikit-learn/scikit-learn/pull/27649
    with ignore_warnings(category=ConvergenceWarning):
        clf_lbfgs.fit(X, y)

    for solver in set(SOLVERS) - set(["lbfgs", "liblinear", "newton-cholesky"]):
        clf = LogisticRegressionCV(solver=solver, **params)
        if solver in ("sag", "saga"):
            clf.set_params(
                tol=1e-18, max_iter=10000, random_state=global_random_seed + 1
            )
        clf.fit(X, y)

        assert_allclose(
            clf.coef_, clf_lbfgs.coef_, rtol=1e-3, err_msg=f"{solver} vs lbfgs"
        )


@pytest.mark.parametrize("problem", ("single", "cv"))
@pytest.mark.parametrize(
    "solver", ("lbfgs", "liblinear", "newton-cg", "newton-cholesky", "sag", "saga")
)
def test_logistic_regression_sample_weights(problem, solver, global_random_seed):
    n_samples_per_cv_group = 200
    n_cv_groups = 3

    X, y = make_classification(
        n_samples=n_samples_per_cv_group * n_cv_groups,
        n_features=5,
        n_informative=3,
        n_classes=2,
        n_redundant=0,
        random_state=global_random_seed,
    )
    rng = np.random.RandomState(global_random_seed)
    sw = np.ones(y.shape[0])

    kw_weighted = {
        "random_state": global_random_seed,
        "fit_intercept": False,
        "max_iter": 100_000 if solver.startswith("sag") else 1_000,
        "tol": 1e-8,
    }
    kw_repeated = kw_weighted.copy()
    sw[:n_samples_per_cv_group] = rng.randint(0, 5, size=n_samples_per_cv_group)
    X_repeated = np.repeat(X, sw.astype(int), axis=0)
    y_repeated = np.repeat(y, sw.astype(int), axis=0)

    if problem == "single":
        LR = LogisticRegression
    elif problem == "cv":
        LR = LogisticRegressionCV
        # We weight the first fold 2 times more.
        groups_weighted = np.concatenate(
            [
                np.full(n_samples_per_cv_group, 0),
                np.full(n_samples_per_cv_group, 1),
                np.full(n_samples_per_cv_group, 2),
            ]
        )
        splits_weighted = list(LeaveOneGroupOut().split(X, groups=groups_weighted))
        kw_weighted.update({"Cs": 100, "cv": splits_weighted})

        groups_repeated = np.repeat(groups_weighted, sw.astype(int), axis=0)
        splits_repeated = list(
            LeaveOneGroupOut().split(X_repeated, groups=groups_repeated)
        )
        kw_repeated.update({"Cs": 100, "cv": splits_repeated})

    clf_sw_weighted = LR(solver=solver, **kw_weighted)
    clf_sw_repeated = LR(solver=solver, **kw_repeated)

    if solver == "lbfgs":
        # lbfgs has convergence issues on the data but this should not impact
        # the quality of the results.
        with warnings.catch_warnings():
            warnings.simplefilter("ignore", ConvergenceWarning)
            clf_sw_weighted.fit(X, y, sample_weight=sw)
            clf_sw_repeated.fit(X_repeated, y_repeated)

    else:
        clf_sw_weighted.fit(X, y, sample_weight=sw)
        clf_sw_repeated.fit(X_repeated, y_repeated)

    if problem == "cv":
        assert_allclose(clf_sw_weighted.scores_[1], clf_sw_repeated.scores_[1])
    assert_allclose(clf_sw_weighted.coef_, clf_sw_repeated.coef_, atol=1e-5)


@pytest.mark.parametrize(
    "solver", ("lbfgs", "newton-cg", "newton-cholesky", "sag", "saga")
)
def test_logistic_regression_solver_class_weights(solver, global_random_seed):
    # Test that passing class_weight as [1, 2] is the same as
    # passing class weight = [1,1] but adjusting sample weights
    # to be 2 for all instances of class 1.

    X, y = make_classification(
        n_samples=300,
        n_features=5,
        n_informative=3,
        n_classes=2,
        random_state=global_random_seed,
    )

    sample_weight = y + 1

    kw_weighted = {
        "random_state": global_random_seed,
        "fit_intercept": False,
        "max_iter": 100_000,
        "tol": 1e-8,
    }
    clf_cw_12 = LogisticRegression(
        solver=solver, class_weight={0: 1, 1: 2}, **kw_weighted
    )
    clf_cw_12.fit(X, y)
    clf_sw_12 = LogisticRegression(solver=solver, **kw_weighted)
    clf_sw_12.fit(X, y, sample_weight=sample_weight)
    assert_allclose(clf_cw_12.coef_, clf_sw_12.coef_, atol=1e-6)


def test_sample_and_class_weight_equivalence_liblinear(global_random_seed):
    # Test the above for l1 penalty and l2 penalty with dual=True.
    # since the patched liblinear code is different.

    X, y = make_classification(
        n_samples=300,
        n_features=5,
        n_informative=3,
        n_classes=2,
        random_state=global_random_seed,
    )

    sample_weight = y + 1

    clf_cw = LogisticRegression(
        solver="liblinear",
        fit_intercept=False,
        class_weight={0: 1, 1: 2},
        penalty="l1",
        max_iter=10_000,
        tol=1e-12,
        random_state=global_random_seed,
    )
    clf_cw.fit(X, y)
    clf_sw = LogisticRegression(
        solver="liblinear",
        fit_intercept=False,
        penalty="l1",
        max_iter=10_000,
        tol=1e-12,
        random_state=global_random_seed,
    )
    clf_sw.fit(X, y, sample_weight)
    assert_allclose(clf_cw.coef_, clf_sw.coef_, atol=1e-10)

    clf_cw = LogisticRegression(
        solver="liblinear",
        fit_intercept=False,
        class_weight={0: 1, 1: 2},
        penalty="l2",
        max_iter=10_000,
        tol=1e-12,
        dual=True,
        random_state=global_random_seed,
    )
    clf_cw.fit(X, y)
    clf_sw = LogisticRegression(
        solver="liblinear",
        fit_intercept=False,
        penalty="l2",
        max_iter=10_000,
        tol=1e-12,
        dual=True,
        random_state=global_random_seed,
    )
    clf_sw.fit(X, y, sample_weight)
    assert_allclose(clf_cw.coef_, clf_sw.coef_, atol=1e-10)


def _compute_class_weight_dictionary(y):
    # helper for returning a dictionary instead of an array
    classes = np.unique(y)
    class_weight = compute_class_weight("balanced", classes=classes, y=y)
    class_weight_dict = dict(zip(classes, class_weight))
    return class_weight_dict


@pytest.mark.parametrize("csr_container", [lambda x: x] + CSR_CONTAINERS)
def test_logistic_regression_class_weights(csr_container):
    # Scale data to avoid convergence warnings with the lbfgs solver
    X_iris = scale(iris.data)
    # Multinomial case: remove 90% of class 0
    X = X_iris[45:, :]
    X = csr_container(X)
    y = iris.target[45:]
    class_weight_dict = _compute_class_weight_dictionary(y)

    for solver in set(SOLVERS) - set(["liblinear", "newton-cholesky"]):
        params = dict(solver=solver, max_iter=1000)
        clf1 = LogisticRegression(class_weight="balanced", **params)
        clf2 = LogisticRegression(class_weight=class_weight_dict, **params)
        clf1.fit(X, y)
        clf2.fit(X, y)
        assert len(clf1.classes_) == 3
        assert_allclose(clf1.coef_, clf2.coef_, rtol=1e-4)
        # Same as appropriate sample_weight.
        sw = np.ones(X.shape[0])
        for c in clf1.classes_:
            sw[y == c] *= class_weight_dict[c]
        clf3 = LogisticRegression(**params).fit(X, y, sample_weight=sw)
        assert_allclose(clf3.coef_, clf2.coef_, rtol=1e-4)

    # Binary case: remove 90% of class 0 and 100% of class 2
    X = X_iris[45:100, :]
    y = iris.target[45:100]
    class_weight_dict = _compute_class_weight_dictionary(y)

    for solver in SOLVERS:
        params = dict(solver=solver, max_iter=1000)
        clf1 = LogisticRegression(class_weight="balanced", **params)
        clf2 = LogisticRegression(class_weight=class_weight_dict, **params)
        clf1.fit(X, y)
        clf2.fit(X, y)
        assert_array_almost_equal(clf1.coef_, clf2.coef_, decimal=6)


def test_logistic_regression_multinomial():
    # Tests for the multinomial option in logistic regression

    # Some basic attributes of Logistic Regression
    n_samples, n_features, n_classes = 50, 20, 3
    X, y = make_classification(
        n_samples=n_samples,
        n_features=n_features,
        n_informative=10,
        n_classes=n_classes,
        random_state=0,
    )

    X = StandardScaler(with_mean=False).fit_transform(X)

    # 'lbfgs' is used as a referenced
    solver = "lbfgs"
    ref_i = LogisticRegression(solver=solver, tol=1e-6)
    ref_w = LogisticRegression(solver=solver, fit_intercept=False, tol=1e-6)
    ref_i.fit(X, y)
    ref_w.fit(X, y)
    assert ref_i.coef_.shape == (n_classes, n_features)
    assert ref_w.coef_.shape == (n_classes, n_features)
    for solver in ["sag", "saga", "newton-cg"]:
        clf_i = LogisticRegression(
            solver=solver,
            random_state=42,
            max_iter=2000,
            tol=1e-7,
        )
        clf_w = LogisticRegression(
            solver=solver,
            random_state=42,
            max_iter=2000,
            tol=1e-7,
            fit_intercept=False,
        )
        clf_i.fit(X, y)
        clf_w.fit(X, y)
        assert clf_i.coef_.shape == (n_classes, n_features)
        assert clf_w.coef_.shape == (n_classes, n_features)

        # Compare solutions between lbfgs and the other solvers
        assert_allclose(ref_i.coef_, clf_i.coef_, rtol=1e-3)
        assert_allclose(ref_w.coef_, clf_w.coef_, rtol=1e-2)
        assert_allclose(ref_i.intercept_, clf_i.intercept_, rtol=1e-3)

    # Test that the path give almost the same results. However since in this
    # case we take the average of the coefs after fitting across all the
    # folds, it need not be exactly the same.
    for solver in ["lbfgs", "newton-cg", "sag", "saga"]:
        clf_path = LogisticRegressionCV(
            solver=solver, max_iter=2000, tol=1e-6, Cs=[1.0]
        )
        clf_path.fit(X, y)
        assert_allclose(clf_path.coef_, ref_i.coef_, rtol=1e-2)
        assert_allclose(clf_path.intercept_, ref_i.intercept_, rtol=1e-2)


def test_liblinear_decision_function_zero():
    # Test negative prediction when decision_function values are zero.
    # Liblinear predicts the positive class when decision_function values
    # are zero. This is a test to verify that we do not do the same.
    # See Issue: https://github.com/scikit-learn/scikit-learn/issues/3600
    # and the PR https://github.com/scikit-learn/scikit-learn/pull/3623
    X, y = make_classification(n_samples=5, n_features=5, random_state=0)
    clf = LogisticRegression(fit_intercept=False, solver="liblinear")
    clf.fit(X, y)

    # Dummy data such that the decision function becomes zero.
    X = np.zeros((5, 5))
    assert_array_equal(clf.predict(X), np.zeros(5))


@pytest.mark.parametrize("csr_container", CSR_CONTAINERS)
def test_liblinear_logregcv_sparse(csr_container):
    # Test LogRegCV with solver='liblinear' works for sparse matrices

    X, y = make_classification(n_samples=10, n_features=5, random_state=0)
    clf = LogisticRegressionCV(solver="liblinear")
    clf.fit(csr_container(X), y)


@pytest.mark.parametrize("csr_container", CSR_CONTAINERS)
def test_saga_sparse(csr_container):
    # Test LogRegCV with solver='liblinear' works for sparse matrices

    X, y = make_classification(n_samples=10, n_features=5, random_state=0)
    clf = LogisticRegressionCV(solver="saga", tol=1e-2)
    clf.fit(csr_container(X), y)


def test_logreg_intercept_scaling_zero():
    # Test that intercept_scaling is ignored when fit_intercept is False

    clf = LogisticRegression(fit_intercept=False)
    clf.fit(X, Y1)
    assert clf.intercept_ == 0.0


def test_logreg_l1():
    # Because liblinear penalizes the intercept and saga does not, we do not
    # fit the intercept to make it possible to compare the coefficients of
    # the two models at convergence.
    rng = np.random.RandomState(42)
    n_samples = 50
    X, y = make_classification(n_samples=n_samples, n_features=20, random_state=0)
    X_noise = rng.normal(size=(n_samples, 3))
    X_constant = np.ones(shape=(n_samples, 2))
    X = np.concatenate((X, X_noise, X_constant), axis=1)
    lr_liblinear = LogisticRegression(
        penalty="l1",
        C=1.0,
        solver="liblinear",
        fit_intercept=False,
        tol=1e-10,
    )
    lr_liblinear.fit(X, y)

    lr_saga = LogisticRegression(
        penalty="l1",
        C=1.0,
        solver="saga",
        fit_intercept=False,
        max_iter=1000,
        tol=1e-10,
    )
    lr_saga.fit(X, y)
    assert_array_almost_equal(lr_saga.coef_, lr_liblinear.coef_)

    # Noise and constant features should be regularized to zero by the l1
    # penalty
    assert_array_almost_equal(lr_liblinear.coef_[0, -5:], np.zeros(5))
    assert_array_almost_equal(lr_saga.coef_[0, -5:], np.zeros(5))


@pytest.mark.parametrize("csr_container", CSR_CONTAINERS)
def test_logreg_l1_sparse_data(csr_container):
    # Because liblinear penalizes the intercept and saga does not, we do not
    # fit the intercept to make it possible to compare the coefficients of
    # the two models at convergence.
    rng = np.random.RandomState(42)
    n_samples = 50
    X, y = make_classification(n_samples=n_samples, n_features=20, random_state=0)
    X_noise = rng.normal(scale=0.1, size=(n_samples, 3))
    X_constant = np.zeros(shape=(n_samples, 2))
    X = np.concatenate((X, X_noise, X_constant), axis=1)
    X[X < 1] = 0
    X = csr_container(X)

    lr_liblinear = LogisticRegression(
        penalty="l1",
        C=1.0,
        solver="liblinear",
        fit_intercept=False,
        tol=1e-10,
    )
    lr_liblinear.fit(X, y)

    lr_saga = LogisticRegression(
        penalty="l1",
        C=1.0,
        solver="saga",
        fit_intercept=False,
        max_iter=1000,
        tol=1e-10,
    )
    lr_saga.fit(X, y)
    assert_array_almost_equal(lr_saga.coef_, lr_liblinear.coef_)
    # Noise and constant features should be regularized to zero by the l1
    # penalty
    assert_array_almost_equal(lr_liblinear.coef_[0, -5:], np.zeros(5))
    assert_array_almost_equal(lr_saga.coef_[0, -5:], np.zeros(5))

    # Check that solving on the sparse and dense data yield the same results
    lr_saga_dense = LogisticRegression(
        penalty="l1",
        C=1.0,
        solver="saga",
        fit_intercept=False,
        max_iter=1000,
        tol=1e-10,
    )
    lr_saga_dense.fit(X.toarray(), y)
    assert_array_almost_equal(lr_saga.coef_, lr_saga_dense.coef_)


@pytest.mark.parametrize("random_seed", [42])
@pytest.mark.parametrize("penalty", ["l1", "l2"])
def test_logistic_regression_cv_refit(random_seed, penalty):
    # Test that when refit=True, logistic regression cv with the saga solver
    # converges to the same solution as logistic regression with a fixed
    # regularization parameter.
    # Internally the LogisticRegressionCV model uses a warm start to refit on
    # the full data model with the optimal C found by CV. As the penalized
    # logistic regression loss is convex, we should still recover exactly
    # the same solution as long as the stopping criterion is strict enough (and
    # that there are no exactly duplicated features when penalty='l1').
    X, y = make_classification(n_samples=100, n_features=20, random_state=random_seed)
    common_params = dict(
        solver="saga",
        penalty=penalty,
        random_state=random_seed,
        max_iter=1000,
        tol=1e-12,
    )
    lr_cv = LogisticRegressionCV(Cs=[1.0], refit=True, **common_params)
    lr_cv.fit(X, y)
    lr = LogisticRegression(C=1.0, **common_params)
    lr.fit(X, y)
    assert_array_almost_equal(lr_cv.coef_, lr.coef_)


def test_logreg_predict_proba_multinomial():
    X, y = make_classification(
        n_samples=10, n_features=20, random_state=0, n_classes=3, n_informative=10
    )

    # Predicted probabilities using the true-entropy loss should give a
    # smaller loss than those using the ovr method.
    clf_multi = LogisticRegression(solver="lbfgs")
    clf_multi.fit(X, y)
    clf_multi_loss = log_loss(y, clf_multi.predict_proba(X))
    clf_ovr = OneVsRestClassifier(LogisticRegression(solver="lbfgs"))
    clf_ovr.fit(X, y)
    clf_ovr_loss = log_loss(y, clf_ovr.predict_proba(X))
    assert clf_ovr_loss > clf_multi_loss

    # Predicted probabilities using the soft-max function should give a
    # smaller loss than those using the logistic function.
    clf_multi_loss = log_loss(y, clf_multi.predict_proba(X))
    clf_wrong_loss = log_loss(y, clf_multi._predict_proba_lr(X))
    assert clf_wrong_loss > clf_multi_loss


# TODO(1.8): remove filterwarnings after the deprecation of multi_class
@pytest.mark.filterwarnings("ignore:.*'multi_class' was deprecated.*:FutureWarning")
@pytest.mark.parametrize("max_iter", np.arange(1, 5))
@pytest.mark.parametrize("multi_class", ["ovr", "multinomial"])
@pytest.mark.parametrize(
    "solver, message",
    [
        (
            "newton-cg",
            "newton-cg failed to converge.* Increase the number of iterations.",
        ),
        (
            "liblinear",
            "Liblinear failed to converge, increase the number of iterations.",
        ),
        ("sag", "The max_iter was reached which means the coef_ did not converge"),
        ("saga", "The max_iter was reached which means the coef_ did not converge"),
        ("lbfgs", "lbfgs failed to converge"),
        ("newton-cholesky", "Newton solver did not converge after [0-9]* iterations"),
    ],
)
def test_max_iter(max_iter, multi_class, solver, message):
    # Test that the maximum number of iteration is reached
    X, y_bin = iris.data, iris.target.copy()
    y_bin[y_bin == 2] = 0

    if solver in ("liblinear",) and multi_class == "multinomial":
        pytest.skip("'multinomial' is not supported by liblinear")
    if solver == "newton-cholesky" and max_iter > 1:
        pytest.skip("solver newton-cholesky might converge very fast")

    lr = LogisticRegression(
        max_iter=max_iter,
        tol=1e-15,
        multi_class=multi_class,
        random_state=0,
        solver=solver,
    )
    with pytest.warns(ConvergenceWarning, match=message):
        lr.fit(X, y_bin)

    assert lr.n_iter_[0] == max_iter


# TODO(1.8): remove filterwarnings after the deprecation of multi_class
@pytest.mark.filterwarnings("ignore:.*'multi_class' was deprecated.*:FutureWarning")
@pytest.mark.filterwarnings(
    "ignore:.*'liblinear' solver for multiclass classification is deprecated.*"
)
@pytest.mark.parametrize("solver", SOLVERS)
def test_n_iter(solver):
    # Test that self.n_iter_ has the correct format.
    X, y = iris.data, iris.target
    if solver == "lbfgs":
        # lbfgs requires scaling to avoid convergence warnings
        X = scale(X)

    n_classes = np.unique(y).shape[0]
    assert n_classes == 3

    # Also generate a binary classification sub-problem.
    y_bin = y.copy()
    y_bin[y_bin == 2] = 0

    n_Cs = 4
    n_cv_fold = 2

    # Binary classification case
    clf = LogisticRegression(tol=1e-2, C=1.0, solver=solver, random_state=42)
    clf.fit(X, y_bin)
    assert clf.n_iter_.shape == (1,)

    clf_cv = LogisticRegressionCV(
        tol=1e-2, solver=solver, Cs=n_Cs, cv=n_cv_fold, random_state=42
    )
    clf_cv.fit(X, y_bin)
    assert clf_cv.n_iter_.shape == (1, n_cv_fold, n_Cs)

    # OvR case
    clf.set_params(multi_class="ovr").fit(X, y)
    assert clf.n_iter_.shape == (n_classes,)

    clf_cv.set_params(multi_class="ovr").fit(X, y)
    assert clf_cv.n_iter_.shape == (n_classes, n_cv_fold, n_Cs)

    # multinomial case
    if solver in ("liblinear",):
        # This solver only supports one-vs-rest multiclass classification.
        return

    # When using the multinomial objective function, there is a single
    # optimization problem to solve for all classes at once:
    clf.set_params(multi_class="multinomial").fit(X, y)
    assert clf.n_iter_.shape == (1,)

    clf_cv.set_params(multi_class="multinomial").fit(X, y)
    assert clf_cv.n_iter_.shape == (1, n_cv_fold, n_Cs)


@pytest.mark.parametrize(
    "solver", sorted(set(SOLVERS) - set(["liblinear", "newton-cholesky"]))
)
@pytest.mark.parametrize("warm_start", (True, False))
@pytest.mark.parametrize("fit_intercept", (True, False))
def test_warm_start(solver, warm_start, fit_intercept):
    # A 1-iteration second fit on same data should give almost same result
    # with warm starting, and quite different result without warm starting.
    # Warm starting does not work with liblinear solver.
    X, y = iris.data, iris.target

    clf = LogisticRegression(
        tol=1e-4,
        warm_start=warm_start,
        solver=solver,
        random_state=42,
        fit_intercept=fit_intercept,
    )
    with ignore_warnings(category=ConvergenceWarning):
        clf.fit(X, y)
        coef_1 = clf.coef_

        clf.max_iter = 1
        clf.fit(X, y)
    cum_diff = np.sum(np.abs(coef_1 - clf.coef_))
    msg = (
        f"Warm starting issue with solver {solver}"
        f"with {fit_intercept=} and {warm_start=}"
    )
    if warm_start:
        assert 2.0 > cum_diff, msg
    else:
        assert cum_diff > 2.0, msg


@pytest.mark.parametrize("csr_container", CSR_CONTAINERS)
def test_saga_vs_liblinear(csr_container):
    iris = load_iris()
    X, y = iris.data, iris.target
    X = np.concatenate([X] * 3)
    y = np.concatenate([y] * 3)

    X_bin = X[y <= 1]
    y_bin = y[y <= 1] * 2 - 1

    X_sparse, y_sparse = make_classification(
        n_samples=50, n_features=20, random_state=0
    )
    X_sparse = csr_container(X_sparse)

    for X, y in ((X_bin, y_bin), (X_sparse, y_sparse)):
        for penalty in ["l1", "l2"]:
            n_samples = X.shape[0]
            # alpha=1e-3 is time consuming
            for alpha in np.logspace(-1, 1, 3):
                saga = LogisticRegression(
                    C=1.0 / (n_samples * alpha),
                    solver="saga",
                    max_iter=200,
                    fit_intercept=False,
                    penalty=penalty,
                    random_state=0,
                    tol=1e-6,
                )

                liblinear = LogisticRegression(
                    C=1.0 / (n_samples * alpha),
                    solver="liblinear",
                    max_iter=200,
                    fit_intercept=False,
                    penalty=penalty,
                    random_state=0,
                    tol=1e-6,
                )

                saga.fit(X, y)
                liblinear.fit(X, y)
                # Convergence for alpha=1e-3 is very slow
                assert_array_almost_equal(saga.coef_, liblinear.coef_, 3)


# TODO(1.8): remove filterwarnings after the deprecation of multi_class
@pytest.mark.filterwarnings("ignore:.*'multi_class' was deprecated.*:FutureWarning")
@pytest.mark.parametrize("multi_class", ["ovr", "multinomial"])
@pytest.mark.parametrize(
    "solver", ["liblinear", "newton-cg", "newton-cholesky", "saga"]
)
@pytest.mark.parametrize("fit_intercept", [False, True])
@pytest.mark.parametrize("csr_container", CSR_CONTAINERS)
def test_dtype_match(solver, multi_class, fit_intercept, csr_container):
    # Test that np.float32 input data is not cast to np.float64 when possible
    # and that the output is approximately the same no matter the input format.

    if solver == "liblinear" and multi_class == "multinomial":
        pytest.skip(f"Solver={solver} does not support multinomial logistic.")

    out32_type = np.float64 if solver == "liblinear" else np.float32

    X_32 = np.array(X).astype(np.float32)
    y_32 = np.array(Y1).astype(np.float32)
    X_64 = np.array(X).astype(np.float64)
    y_64 = np.array(Y1).astype(np.float64)
    X_sparse_32 = csr_container(X, dtype=np.float32)
    X_sparse_64 = csr_container(X, dtype=np.float64)
    solver_tol = 5e-4

    lr_templ = LogisticRegression(
        solver=solver,
        multi_class=multi_class,
        random_state=42,
        tol=solver_tol,
        fit_intercept=fit_intercept,
    )

    # Check 32-bit type consistency
    lr_32 = clone(lr_templ)
    lr_32.fit(X_32, y_32)
    assert lr_32.coef_.dtype == out32_type

    # Check 32-bit type consistency with sparsity
    lr_32_sparse = clone(lr_templ)
    lr_32_sparse.fit(X_sparse_32, y_32)
    assert lr_32_sparse.coef_.dtype == out32_type

    # Check 64-bit type consistency
    lr_64 = clone(lr_templ)
    lr_64.fit(X_64, y_64)
    assert lr_64.coef_.dtype == np.float64

    # Check 64-bit type consistency with sparsity
    lr_64_sparse = clone(lr_templ)
    lr_64_sparse.fit(X_sparse_64, y_64)
    assert lr_64_sparse.coef_.dtype == np.float64

    # solver_tol bounds the norm of the loss gradient
    # dw ~= inv(H)*grad ==> |dw| ~= |inv(H)| * solver_tol, where H - hessian
    #
    # See https://github.com/scikit-learn/scikit-learn/pull/13645
    #
    # with  Z = np.hstack((np.ones((3,1)), np.array(X)))
    # In [8]: np.linalg.norm(np.diag([0,2,2]) + np.linalg.inv((Z.T @ Z)/4))
    # Out[8]: 1.7193336918135917

    # factor of 2 to get the ball diameter
    atol = 2 * 1.72 * solver_tol
    if os.name == "nt" and _IS_32BIT:
        # FIXME
        atol = 1e-2

    # Check accuracy consistency
    assert_allclose(lr_32.coef_, lr_64.coef_.astype(np.float32), atol=atol)

    if solver == "saga" and fit_intercept:
        # FIXME: SAGA on sparse data fits the intercept inaccurately with the
        # default tol and max_iter parameters.
        atol = 1e-1

    assert_allclose(lr_32.coef_, lr_32_sparse.coef_, atol=atol)
    assert_allclose(lr_64.coef_, lr_64_sparse.coef_, atol=atol)


def test_warm_start_converge_LR():
    # Test to see that the logistic regression converges on warm start,
    # with multi_class='multinomial'. Non-regressive test for #10836

    rng = np.random.RandomState(0)
    X = np.concatenate((rng.randn(100, 2) + [1, 1], rng.randn(100, 2)))
    y = np.array([1] * 100 + [-1] * 100)
    lr_no_ws = LogisticRegression(solver="sag", warm_start=False, random_state=0)
    lr_ws = LogisticRegression(solver="sag", warm_start=True, random_state=0)

    lr_no_ws_loss = log_loss(y, lr_no_ws.fit(X, y).predict_proba(X))
    for i in range(5):
        lr_ws.fit(X, y)
    lr_ws_loss = log_loss(y, lr_ws.predict_proba(X))
    assert_allclose(lr_no_ws_loss, lr_ws_loss, rtol=1e-5)


def test_elastic_net_coeffs():
    # make sure elasticnet penalty gives different coefficients from l1 and l2
    # with saga solver (l1_ratio different from 0 or 1)
    X, y = make_classification(random_state=0)

    C = 2.0
    l1_ratio = 0.5
    coeffs = list()
    for penalty, ratio in (("elasticnet", l1_ratio), ("l1", None), ("l2", None)):
        lr = LogisticRegression(
            penalty=penalty,
            C=C,
            solver="saga",
            random_state=0,
            l1_ratio=ratio,
            tol=1e-3,
            max_iter=200,
        )
        lr.fit(X, y)
        coeffs.append(lr.coef_)

    elastic_net_coeffs, l1_coeffs, l2_coeffs = coeffs
    # make sure coeffs differ by at least .1
    assert not np.allclose(elastic_net_coeffs, l1_coeffs, rtol=0, atol=0.1)
    assert not np.allclose(elastic_net_coeffs, l2_coeffs, rtol=0, atol=0.1)
    assert not np.allclose(l2_coeffs, l1_coeffs, rtol=0, atol=0.1)


@pytest.mark.parametrize("C", [0.001, 0.1, 1, 10, 100, 1000, 1e6])
@pytest.mark.parametrize("penalty, l1_ratio", [("l1", 1), ("l2", 0)])
def test_elastic_net_l1_l2_equivalence(C, penalty, l1_ratio):
    # Make sure elasticnet is equivalent to l1 when l1_ratio=1 and to l2 when
    # l1_ratio=0.
    X, y = make_classification(random_state=0)

    lr_enet = LogisticRegression(
        penalty="elasticnet",
        C=C,
        l1_ratio=l1_ratio,
        solver="saga",
        random_state=0,
        tol=1e-2,
    )
    lr_expected = LogisticRegression(
        penalty=penalty, C=C, solver="saga", random_state=0, tol=1e-2
    )
    lr_enet.fit(X, y)
    lr_expected.fit(X, y)

    assert_array_almost_equal(lr_enet.coef_, lr_expected.coef_)


@pytest.mark.parametrize("C", [0.001, 1, 100, 1e6])
def test_elastic_net_vs_l1_l2(C):
    # Make sure that elasticnet with grid search on l1_ratio gives same or
    # better results than just l1 or just l2.

    X, y = make_classification(500, random_state=0)
    X_train, X_test, y_train, y_test = train_test_split(X, y, random_state=0)

    param_grid = {"l1_ratio": np.linspace(0, 1, 5)}

    enet_clf = LogisticRegression(
        penalty="elasticnet", C=C, solver="saga", random_state=0, tol=1e-2
    )
    gs = GridSearchCV(enet_clf, param_grid, refit=True)

    l1_clf = LogisticRegression(
        penalty="l1", C=C, solver="saga", random_state=0, tol=1e-2
    )
    l2_clf = LogisticRegression(
        penalty="l2", C=C, solver="saga", random_state=0, tol=1e-2
    )

    for clf in (gs, l1_clf, l2_clf):
        clf.fit(X_train, y_train)

    assert gs.score(X_test, y_test) >= l1_clf.score(X_test, y_test)
    assert gs.score(X_test, y_test) >= l2_clf.score(X_test, y_test)


@pytest.mark.parametrize("C", np.logspace(-3, 2, 4))
@pytest.mark.parametrize("l1_ratio", [0.1, 0.5, 0.9])
def test_LogisticRegression_elastic_net_objective(C, l1_ratio):
    # Check that training with a penalty matching the objective leads
    # to a lower objective.
    # Here we train a logistic regression with l2 (a) and elasticnet (b)
    # penalties, and compute the elasticnet objective. That of a should be
    # greater than that of b (both objectives are convex).
    X, y = make_classification(
        n_samples=1000,
        n_classes=2,
        n_features=20,
        n_informative=10,
        n_redundant=0,
        n_repeated=0,
        random_state=0,
    )
    X = scale(X)

    lr_enet = LogisticRegression(
        penalty="elasticnet",
        solver="saga",
        random_state=0,
        C=C,
        l1_ratio=l1_ratio,
        fit_intercept=False,
    )
    lr_l2 = LogisticRegression(
        penalty="l2", solver="saga", random_state=0, C=C, fit_intercept=False
    )
    lr_enet.fit(X, y)
    lr_l2.fit(X, y)

    def enet_objective(lr):
        coef = lr.coef_.ravel()
        obj = C * log_loss(y, lr.predict_proba(X))
        obj += l1_ratio * np.sum(np.abs(coef))
        obj += (1.0 - l1_ratio) * 0.5 * np.dot(coef, coef)
        return obj

    assert enet_objective(lr_enet) < enet_objective(lr_l2)


@pytest.mark.parametrize("n_classes", (2, 3))
def test_LogisticRegressionCV_GridSearchCV_elastic_net(n_classes):
    # make sure LogisticRegressionCV gives same best params (l1 and C) as
    # GridSearchCV when penalty is elasticnet

    X, y = make_classification(
        n_samples=100, n_classes=n_classes, n_informative=3, random_state=0
    )

    cv = StratifiedKFold(5)

    l1_ratios = np.linspace(0, 1, 3)
    Cs = np.logspace(-4, 4, 3)

    lrcv = LogisticRegressionCV(
        penalty="elasticnet",
        Cs=Cs,
        solver="saga",
        cv=cv,
        l1_ratios=l1_ratios,
        random_state=0,
        tol=1e-2,
    )
    lrcv.fit(X, y)

    param_grid = {"C": Cs, "l1_ratio": l1_ratios}
    lr = LogisticRegression(
        penalty="elasticnet",
        solver="saga",
        random_state=0,
        tol=1e-2,
    )
    gs = GridSearchCV(lr, param_grid, cv=cv)
    gs.fit(X, y)

    assert gs.best_params_["l1_ratio"] == lrcv.l1_ratio_[0]
    assert gs.best_params_["C"] == lrcv.C_[0]


# TODO(1.8): remove filterwarnings after the deprecation of multi_class
# Maybe remove whole test after removal of the deprecated multi_class.
@pytest.mark.filterwarnings("ignore:.*'multi_class' was deprecated.*:FutureWarning")
def test_LogisticRegressionCV_GridSearchCV_elastic_net_ovr():
    # make sure LogisticRegressionCV gives same best params (l1 and C) as
    # GridSearchCV when penalty is elasticnet and multiclass is ovr. We can't
    # compare best_params like in the previous test because
    # LogisticRegressionCV with multi_class='ovr' will have one C and one
    # l1_param for each class, while LogisticRegression will share the
    # parameters over the *n_classes* classifiers.

    X, y = make_classification(
        n_samples=100, n_classes=3, n_informative=3, random_state=0
    )
    X_train, X_test, y_train, y_test = train_test_split(X, y, random_state=0)
    cv = StratifiedKFold(5)

    l1_ratios = np.linspace(0, 1, 3)
    Cs = np.logspace(-4, 4, 3)

    lrcv = LogisticRegressionCV(
        penalty="elasticnet",
        Cs=Cs,
        solver="saga",
        cv=cv,
        l1_ratios=l1_ratios,
        random_state=0,
        multi_class="ovr",
        tol=1e-2,
    )
    lrcv.fit(X_train, y_train)

    param_grid = {"C": Cs, "l1_ratio": l1_ratios}
    lr = LogisticRegression(
        penalty="elasticnet",
        solver="saga",
        random_state=0,
        multi_class="ovr",
        tol=1e-2,
    )
    gs = GridSearchCV(lr, param_grid, cv=cv)
    gs.fit(X_train, y_train)

    # Check that predictions are 80% the same
    assert (lrcv.predict(X_train) == gs.predict(X_train)).mean() >= 0.8
    assert (lrcv.predict(X_test) == gs.predict(X_test)).mean() >= 0.8


# TODO(1.8): remove filterwarnings after the deprecation of multi_class
@pytest.mark.filterwarnings("ignore:.*'multi_class' was deprecated.*:FutureWarning")
@pytest.mark.parametrize("penalty", ("l2", "elasticnet"))
@pytest.mark.parametrize("multi_class", ("ovr", "multinomial", "auto"))
def test_LogisticRegressionCV_no_refit(penalty, multi_class):
    # Test LogisticRegressionCV attribute shapes when refit is False

    n_classes = 3
    n_features = 20
    X, y = make_classification(
        n_samples=200,
        n_classes=n_classes,
        n_informative=n_classes,
        n_features=n_features,
        random_state=0,
    )

    Cs = np.logspace(-4, 4, 3)
    if penalty == "elasticnet":
        l1_ratios = np.linspace(0, 1, 2)
    else:
        l1_ratios = None

    lrcv = LogisticRegressionCV(
        penalty=penalty,
        Cs=Cs,
        solver="saga",
        l1_ratios=l1_ratios,
        random_state=0,
        multi_class=multi_class,
        tol=1e-2,
        refit=False,
    )
    lrcv.fit(X, y)
    assert lrcv.C_.shape == (n_classes,)
    assert lrcv.l1_ratio_.shape == (n_classes,)
    assert lrcv.coef_.shape == (n_classes, n_features)


# TODO(1.8): remove filterwarnings after the deprecation of multi_class
# Remove multi_class an change first element of the expected n_iter_.shape from
# n_classes to 1 (according to the docstring).
@pytest.mark.filterwarnings("ignore:.*'multi_class' was deprecated.*:FutureWarning")
def test_LogisticRegressionCV_elasticnet_attribute_shapes():
    # Make sure the shapes of scores_ and coefs_paths_ attributes are correct
    # when using elasticnet (added one dimension for l1_ratios)

    n_classes = 3
    n_features = 20
    X, y = make_classification(
        n_samples=200,
        n_classes=n_classes,
        n_informative=n_classes,
        n_features=n_features,
        random_state=0,
    )

    Cs = np.logspace(-4, 4, 3)
    l1_ratios = np.linspace(0, 1, 2)

    n_folds = 2
    lrcv = LogisticRegressionCV(
        penalty="elasticnet",
        Cs=Cs,
        solver="saga",
        cv=n_folds,
        l1_ratios=l1_ratios,
        multi_class="ovr",
        random_state=0,
        tol=1e-2,
    )
    lrcv.fit(X, y)
    coefs_paths = np.asarray(list(lrcv.coefs_paths_.values()))
    assert coefs_paths.shape == (
        n_classes,
        n_folds,
        Cs.size,
        l1_ratios.size,
        n_features + 1,
    )
    scores = np.asarray(list(lrcv.scores_.values()))
    assert scores.shape == (n_classes, n_folds, Cs.size, l1_ratios.size)

    assert lrcv.n_iter_.shape == (n_classes, n_folds, Cs.size, l1_ratios.size)


def test_l1_ratio_non_elasticnet():
    msg = (
        r"l1_ratio parameter is only used when penalty is"
        r" 'elasticnet'\. Got \(penalty=l1\)"
    )
    with pytest.warns(UserWarning, match=msg):
        LogisticRegression(penalty="l1", solver="saga", l1_ratio=0.5).fit(X, Y1)


@pytest.mark.parametrize("C", np.logspace(-3, 2, 4))
@pytest.mark.parametrize("l1_ratio", [0.1, 0.5, 0.9])
def test_elastic_net_versus_sgd(C, l1_ratio):
    # Compare elasticnet penalty in LogisticRegression() and SGD(loss='log')
    n_samples = 500
    X, y = make_classification(
        n_samples=n_samples,
        n_classes=2,
        n_features=5,
        n_informative=5,
        n_redundant=0,
        n_repeated=0,
        random_state=1,
    )
    X = scale(X)

    sgd = SGDClassifier(
        penalty="elasticnet",
        random_state=1,
        fit_intercept=False,
        tol=None,
        max_iter=2000,
        l1_ratio=l1_ratio,
        alpha=1.0 / C / n_samples,
        loss="log_loss",
    )
    log = LogisticRegression(
        penalty="elasticnet",
        random_state=1,
        fit_intercept=False,
        tol=1e-5,
        max_iter=1000,
        l1_ratio=l1_ratio,
        C=C,
        solver="saga",
    )

    sgd.fit(X, y)
    log.fit(X, y)
    assert_array_almost_equal(sgd.coef_, log.coef_, decimal=1)


def test_logistic_regression_path_coefs_multinomial():
    # Make sure that the returned coefs by logistic_regression_path when
    # multi_class='multinomial' don't override each other (used to be a
    # bug).
    X, y = make_classification(
        n_samples=200,
        n_classes=3,
        n_informative=2,
        n_redundant=0,
        n_clusters_per_class=1,
        random_state=0,
        n_features=2,
    )
    Cs = [0.00001, 1, 10000]
    coefs, _, _ = _logistic_regression_path(
        X,
        y,
        penalty="l1",
        Cs=Cs,
        solver="saga",
        random_state=0,
        multi_class="multinomial",
    )

    with pytest.raises(AssertionError):
        assert_array_almost_equal(coefs[0], coefs[1], decimal=1)
    with pytest.raises(AssertionError):
        assert_array_almost_equal(coefs[0], coefs[2], decimal=1)
    with pytest.raises(AssertionError):
        assert_array_almost_equal(coefs[1], coefs[2], decimal=1)


# TODO(1.8): remove filterwarnings after the deprecation of multi_class
@pytest.mark.filterwarnings("ignore:.*'multi_class' was deprecated.*:FutureWarning")
@pytest.mark.filterwarnings(
    "ignore:.*'liblinear' solver for multiclass classification is deprecated.*"
)
@pytest.mark.parametrize(
    "est",
    [
        LogisticRegression(random_state=0, max_iter=500),
        LogisticRegressionCV(random_state=0, cv=3, Cs=3, tol=1e-3, max_iter=500),
    ],
    ids=lambda x: x.__class__.__name__,
)
@pytest.mark.parametrize("solver", SOLVERS)
def test_logistic_regression_multi_class_auto(est, solver):
    # check multi_class='auto' => multi_class='ovr'
    # iff binary y or liblinear

    def fit(X, y, **kw):
        return clone(est).set_params(**kw).fit(X, y)

    scaled_data = scale(iris.data)
    X = scaled_data[::10]
    X2 = scaled_data[1::10]
    y_multi = iris.target[::10]
    y_bin = y_multi == 0
    est_auto_bin = fit(X, y_bin, multi_class="auto", solver=solver)
    est_ovr_bin = fit(X, y_bin, multi_class="ovr", solver=solver)
    assert_allclose(est_auto_bin.coef_, est_ovr_bin.coef_)
    assert_allclose(est_auto_bin.predict_proba(X2), est_ovr_bin.predict_proba(X2))

    est_auto_multi = fit(X, y_multi, multi_class="auto", solver=solver)
    if solver == "liblinear":
        est_ovr_multi = fit(X, y_multi, multi_class="ovr", solver=solver)
        assert_allclose(est_auto_multi.coef_, est_ovr_multi.coef_)
        assert_allclose(
            est_auto_multi.predict_proba(X2), est_ovr_multi.predict_proba(X2)
        )
    else:
        est_multi_multi = fit(X, y_multi, multi_class="multinomial", solver=solver)
        assert_allclose(est_auto_multi.coef_, est_multi_multi.coef_)
        assert_allclose(
            est_auto_multi.predict_proba(X2), est_multi_multi.predict_proba(X2)
        )

        # Make sure multi_class='ovr' is distinct from ='multinomial'
        assert not np.allclose(
            est_auto_bin.coef_,
            fit(X, y_bin, multi_class="multinomial", solver=solver).coef_,
        )
        assert not np.allclose(
            est_auto_bin.coef_,
            fit(X, y_multi, multi_class="multinomial", solver=solver).coef_,
        )


@pytest.mark.parametrize("solver", sorted(set(SOLVERS) - set(["liblinear"])))
def test_penalty_none(solver):
    # - Make sure warning is raised if penalty=None and C is set to a
    #   non-default value.
    # - Make sure setting penalty=None is equivalent to setting C=np.inf with
    #   l2 penalty.
    X, y = make_classification(n_samples=1000, n_redundant=0, random_state=0)

    msg = "Setting penalty=None will ignore the C"
    lr = LogisticRegression(penalty=None, solver=solver, C=4)
    with pytest.warns(UserWarning, match=msg):
        lr.fit(X, y)

    lr_none = LogisticRegression(penalty=None, solver=solver, random_state=0)
    lr_l2_C_inf = LogisticRegression(
        penalty="l2", C=np.inf, solver=solver, random_state=0
    )
    pred_none = lr_none.fit(X, y).predict(X)
    pred_l2_C_inf = lr_l2_C_inf.fit(X, y).predict(X)
    assert_array_equal(pred_none, pred_l2_C_inf)


@pytest.mark.parametrize(
    "params",
    [
        {"penalty": "l1", "dual": False, "tol": 1e-6, "max_iter": 1000},
        {"penalty": "l2", "dual": True, "tol": 1e-12, "max_iter": 1000},
        {"penalty": "l2", "dual": False, "tol": 1e-12, "max_iter": 1000},
    ],
)
def test_logisticregression_liblinear_sample_weight(params):
    # check that we support sample_weight with liblinear in all possible cases:
    # l1-primal, l2-primal, l2-dual
    X = np.array(
        [
            [1, 3],
            [1, 3],
            [1, 3],
            [1, 3],
            [2, 1],
            [2, 1],
            [2, 1],
            [2, 1],
            [3, 3],
            [3, 3],
            [3, 3],
            [3, 3],
            [4, 1],
            [4, 1],
            [4, 1],
            [4, 1],
        ],
        dtype=np.dtype("float"),
    )
    y = np.array(
        [1, 1, 1, 1, 2, 2, 2, 2, 1, 1, 1, 1, 2, 2, 2, 2], dtype=np.dtype("int")
    )

    X2 = np.vstack([X, X])
    y2 = np.hstack([y, 3 - y])
    sample_weight = np.ones(shape=len(y) * 2)
    sample_weight[len(y) :] = 0
    X2, y2, sample_weight = shuffle(X2, y2, sample_weight, random_state=0)

    base_clf = LogisticRegression(solver="liblinear", random_state=42)
    base_clf.set_params(**params)
    clf_no_weight = clone(base_clf).fit(X, y)
    clf_with_weight = clone(base_clf).fit(X2, y2, sample_weight=sample_weight)

    for method in ("predict", "predict_proba", "decision_function"):
        X_clf_no_weight = getattr(clf_no_weight, method)(X)
        X_clf_with_weight = getattr(clf_with_weight, method)(X)
        assert_allclose(X_clf_no_weight, X_clf_with_weight)


def test_scores_attribute_layout_elasticnet():
    # Non regression test for issue #14955.
    # when penalty is elastic net the scores_ attribute has shape
    # (n_classes, n_Cs, n_l1_ratios)
    # We here make sure that the second dimension indeed corresponds to Cs and
    # the third dimension corresponds to l1_ratios.

    X, y = make_classification(n_samples=1000, random_state=0)
    cv = StratifiedKFold(n_splits=5)

    l1_ratios = [0.1, 0.9]
    Cs = [0.1, 1, 10]

    lrcv = LogisticRegressionCV(
        penalty="elasticnet",
        solver="saga",
        l1_ratios=l1_ratios,
        Cs=Cs,
        cv=cv,
        random_state=0,
        max_iter=250,
        tol=1e-3,
    )
    lrcv.fit(X, y)

    avg_scores_lrcv = lrcv.scores_[1].mean(axis=0)  # average over folds

    for i, C in enumerate(Cs):
        for j, l1_ratio in enumerate(l1_ratios):
            lr = LogisticRegression(
                penalty="elasticnet",
                solver="saga",
                C=C,
                l1_ratio=l1_ratio,
                random_state=0,
                max_iter=250,
                tol=1e-3,
            )

            avg_score_lr = cross_val_score(lr, X, y, cv=cv).mean()
            assert avg_scores_lrcv[i, j] == pytest.approx(avg_score_lr)


# TODO(1.8): remove filterwarnings after the deprecation of multi_class
@pytest.mark.filterwarnings("ignore:.*'multi_class' was deprecated.*:FutureWarning")
@pytest.mark.parametrize("solver", ["lbfgs", "newton-cg", "newton-cholesky"])
@pytest.mark.parametrize("fit_intercept", [False, True])
def test_multinomial_identifiability_on_iris(solver, fit_intercept):
    """Test that the multinomial classification is identifiable.

    A multinomial with c classes can be modeled with
    probability_k = exp(X@coef_k) / sum(exp(X@coef_l), l=1..c) for k=1..c.
    This is not identifiable, unless one chooses a further constraint.
    According to [1], the maximum of the L2 penalized likelihood automatically
    satisfies the symmetric constraint:
    sum(coef_k, k=1..c) = 0

    Further details can be found in [2].

    Reference
    ---------
    .. [1] :doi:`Zhu, Ji and Trevor J. Hastie. "Classification of gene microarrays by
           penalized logistic regression". Biostatistics 5 3 (2004): 427-43.
           <10.1093/biostatistics/kxg046>`

    .. [2] :arxiv:`Noah Simon and Jerome Friedman and Trevor Hastie. (2013)
           "A Blockwise Descent Algorithm for Group-penalized Multiresponse and
           Multinomial Regression". <1311.6529>`
    """
    # Test logistic regression with the iris dataset
    n_samples, n_features = iris.data.shape
    target = iris.target_names[iris.target]

    clf = LogisticRegression(
        C=len(iris.data),
        solver="lbfgs",
        fit_intercept=fit_intercept,
    )
    # Scaling X to ease convergence.
    X_scaled = scale(iris.data)
    clf.fit(X_scaled, target)

    # axis=0 is sum over classes
    assert_allclose(clf.coef_.sum(axis=0), 0, atol=1e-10)
    if fit_intercept:
        assert clf.intercept_.sum(axis=0) == pytest.approx(0, abs=1e-11)


# TODO(1.8): remove filterwarnings after the deprecation of multi_class
@pytest.mark.filterwarnings("ignore:.*'multi_class' was deprecated.*:FutureWarning")
@pytest.mark.parametrize("multi_class", ["ovr", "multinomial", "auto"])
@pytest.mark.parametrize("class_weight", [{0: 1.0, 1: 10.0, 2: 1.0}, "balanced"])
def test_sample_weight_not_modified(multi_class, class_weight):
    X, y = load_iris(return_X_y=True)
    n_features = len(X)
    W = np.ones(n_features)
    W[: n_features // 2] = 2

    expected = W.copy()

    clf = LogisticRegression(
        random_state=0, class_weight=class_weight, max_iter=200, multi_class=multi_class
    )
    clf.fit(X, y, sample_weight=W)
    assert_allclose(expected, W)


@pytest.mark.parametrize("solver", SOLVERS)
@pytest.mark.parametrize("csr_container", CSR_CONTAINERS)
def test_large_sparse_matrix(solver, global_random_seed, csr_container):
    # Solvers either accept large sparse matrices, or raise helpful error.
    # Non-regression test for pull-request #21093.

    # generate sparse matrix with int64 indices
    X = csr_container(sparse.rand(20, 10, random_state=global_random_seed))
    for attr in ["indices", "indptr"]:
        setattr(X, attr, getattr(X, attr).astype("int64"))
    rng = np.random.RandomState(global_random_seed)
    y = rng.randint(2, size=X.shape[0])

    if solver in ["liblinear", "sag", "saga"]:
        msg = "Only sparse matrices with 32-bit integer indices"
        with pytest.raises(ValueError, match=msg):
            LogisticRegression(solver=solver).fit(X, y)
    else:
        LogisticRegression(solver=solver).fit(X, y)


def test_single_feature_newton_cg():
    # Test that Newton-CG works with a single feature and intercept.
    # Non-regression test for issue #23605.

    X = np.array([[0.5, 0.65, 1.1, 1.25, 0.8, 0.54, 0.95, 0.7]]).T
    y = np.array([1, 1, 0, 0, 1, 1, 0, 1])
    assert X.shape[1] == 1
    LogisticRegression(solver="newton-cg", fit_intercept=True).fit(X, y)


def test_liblinear_not_stuck():
    # Non-regression https://github.com/scikit-learn/scikit-learn/issues/18264
    X = iris.data.copy()
    y = iris.target.copy()
    X = X[y != 2]
    y = y[y != 2]
    X_prep = StandardScaler().fit_transform(X)

    C = l1_min_c(X, y, loss="log") * 10 ** (10 / 29)
    clf = LogisticRegression(
        penalty="l1",
        solver="liblinear",
        tol=1e-6,
        max_iter=100,
        intercept_scaling=10000.0,
        random_state=0,
        C=C,
    )

    # test that the fit does not raise a ConvergenceWarning
    with warnings.catch_warnings():
        warnings.simplefilter("error", ConvergenceWarning)
        clf.fit(X_prep, y)


@config_context(enable_metadata_routing=True)
def test_lr_cv_scores_differ_when_sample_weight_is_requested():
    """Test that `sample_weight` is correctly passed to the scorer in
    `LogisticRegressionCV.fit` and `LogisticRegressionCV.score` by
    checking the difference in scores with the case when `sample_weight`
    is not requested.
    """
    rng = np.random.RandomState(10)
    X, y = make_classification(n_samples=10, random_state=rng)
    X_t, y_t = make_classification(n_samples=10, random_state=rng)
    sample_weight = np.ones(len(y))
    sample_weight[: len(y) // 2] = 2
    kwargs = {"sample_weight": sample_weight}

    scorer1 = get_scorer("accuracy")
    lr_cv1 = LogisticRegressionCV(scoring=scorer1)
    lr_cv1.fit(X, y, **kwargs)

    scorer2 = get_scorer("accuracy")
    scorer2.set_score_request(sample_weight=True)
    lr_cv2 = LogisticRegressionCV(scoring=scorer2)
    lr_cv2.fit(X, y, **kwargs)

    assert not np.allclose(lr_cv1.scores_[1], lr_cv2.scores_[1])

    score_1 = lr_cv1.score(X_t, y_t, **kwargs)
    score_2 = lr_cv2.score(X_t, y_t, **kwargs)

    assert not np.allclose(score_1, score_2)


def test_lr_cv_scores_without_enabling_metadata_routing():
    """Test that `sample_weight` is passed correctly to the scorer in
    `LogisticRegressionCV.fit` and `LogisticRegressionCV.score` even
    when `enable_metadata_routing=False`
    """
    rng = np.random.RandomState(10)
    X, y = make_classification(n_samples=10, random_state=rng)
    X_t, y_t = make_classification(n_samples=10, random_state=rng)
    sample_weight = np.ones(len(y))
    sample_weight[: len(y) // 2] = 2
    kwargs = {"sample_weight": sample_weight}

    with config_context(enable_metadata_routing=False):
        scorer1 = get_scorer("accuracy")
        lr_cv1 = LogisticRegressionCV(scoring=scorer1)
        lr_cv1.fit(X, y, **kwargs)
        score_1 = lr_cv1.score(X_t, y_t, **kwargs)

    with config_context(enable_metadata_routing=True):
        scorer2 = get_scorer("accuracy")
        scorer2.set_score_request(sample_weight=True)
        lr_cv2 = LogisticRegressionCV(scoring=scorer2)
        lr_cv2.fit(X, y, **kwargs)
        score_2 = lr_cv2.score(X_t, y_t, **kwargs)

    assert_allclose(lr_cv1.scores_[1], lr_cv2.scores_[1])
    assert_allclose(score_1, score_2)


@pytest.mark.parametrize("solver", SOLVERS)
def test_zero_max_iter(solver):
    # Make sure we can inspect the state of LogisticRegression right after
    # initialization (before the first weight update).
    X, y = load_iris(return_X_y=True)
    y = y == 2
    with ignore_warnings(category=ConvergenceWarning):
        clf = LogisticRegression(solver=solver, max_iter=0).fit(X, y)
    if solver not in ["saga", "sag"]:
        # XXX: sag and saga have n_iter_ = [1]...
        assert clf.n_iter_ == 0

    if solver != "lbfgs":
        # XXX: lbfgs has already started to update the coefficients...
        assert_allclose(clf.coef_, np.zeros_like(clf.coef_))
        assert_allclose(
            clf.decision_function(X),
            np.full(shape=X.shape[0], fill_value=clf.intercept_),
        )
        assert_allclose(
            clf.predict_proba(X),
            np.full(shape=(X.shape[0], 2), fill_value=0.5),
        )
    assert clf.score(X, y) < 0.7


def test_passing_params_without_enabling_metadata_routing():
    """Test that the right error message is raised when metadata params
    are passed while not supported when `enable_metadata_routing=False`."""
    X, y = make_classification(n_samples=10, random_state=0)
    lr_cv = LogisticRegressionCV()
    msg = "is only supported if enable_metadata_routing=True"

    with config_context(enable_metadata_routing=False):
        params = {"extra_param": 1.0}

        with pytest.raises(ValueError, match=msg):
            lr_cv.fit(X, y, **params)

        with pytest.raises(ValueError, match=msg):
            lr_cv.score(X, y, **params)


# TODO(1.8): remove
def test_multi_class_deprecated():
    """Check `multi_class` parameter deprecated."""
    X, y = make_classification(n_classes=3, n_samples=50, n_informative=6)
    lr = LogisticRegression(multi_class="ovr")
    msg = "'multi_class' was deprecated"
    with pytest.warns(FutureWarning, match=msg):
        lr.fit(X, y)

    lrCV = LogisticRegressionCV(multi_class="ovr")
    with pytest.warns(FutureWarning, match=msg):
        lrCV.fit(X, y)

    # Special warning for "binary multinomial"
    X, y = make_classification(n_classes=2, n_samples=50, n_informative=6)
    lr = LogisticRegression(multi_class="multinomial")
    msg = "'multi_class' was deprecated.*binary problems"
    with pytest.warns(FutureWarning, match=msg):
        lr.fit(X, y)

    lrCV = LogisticRegressionCV(multi_class="multinomial")
    with pytest.warns(FutureWarning, match=msg):
        lrCV.fit(X, y)


def test_newton_cholesky_fallback_to_lbfgs(global_random_seed):
    # Wide data matrix should lead to a rank-deficient Hessian matrix
    # hence make the Newton-Cholesky solver raise a warning and fallback to
    # lbfgs.
    X, y = make_classification(
        n_samples=10, n_features=20, random_state=global_random_seed
    )
    C = 1e30  # very high C to nearly disable regularization

    # Check that LBFGS can converge without any warning on this problem.
    lr_lbfgs = LogisticRegression(solver="lbfgs", C=C)
    with warnings.catch_warnings():
        warnings.simplefilter("error")
        lr_lbfgs.fit(X, y)
        n_iter_lbfgs = lr_lbfgs.n_iter_[0]

    assert n_iter_lbfgs >= 1

    # Check that the Newton-Cholesky solver raises a warning and falls back to
    # LBFGS. This should converge with the same number of iterations as the
    # above call of lbfgs since the Newton-Cholesky triggers the fallback
    # before completing the first iteration, for the problem setting at hand.
    lr_nc = LogisticRegression(solver="newton-cholesky", C=C)
    with ignore_warnings(category=LinAlgWarning):
        lr_nc.fit(X, y)
        n_iter_nc = lr_nc.n_iter_[0]

    assert n_iter_nc == n_iter_lbfgs

    # Trying to fit the same model again with a small iteration budget should
    # therefore raise a ConvergenceWarning:
    lr_nc_limited = LogisticRegression(
        solver="newton-cholesky", C=C, max_iter=n_iter_lbfgs - 1
    )
    with ignore_warnings(category=LinAlgWarning):
        with pytest.warns(ConvergenceWarning, match="lbfgs failed to converge"):
            lr_nc_limited.fit(X, y)
            n_iter_nc_limited = lr_nc_limited.n_iter_[0]

    assert n_iter_nc_limited == lr_nc_limited.max_iter - 1


# TODO(1.8): check for an error instead
@pytest.mark.parametrize("Estimator", [LogisticRegression, LogisticRegressionCV])
def test_liblinear_multiclass_warning(Estimator):
    """Check that liblinear warns on multiclass problems."""
    msg = (
        "Using the 'liblinear' solver for multiclass classification is "
        "deprecated. An error will be raised in 1.8. Either use another "
        "solver which supports the multinomial loss or wrap the estimator "
        "in a OneVsRestClassifier to keep applying a one-versus-rest "
        "scheme."
    )
    with pytest.warns(FutureWarning, match=msg):
        Estimator(solver="liblinear").fit(iris.data, iris.target)
