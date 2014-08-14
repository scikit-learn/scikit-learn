import numpy as np
import scipy.sparse as sp
from scipy import linalg

from sklearn.utils.testing import assert_true
from sklearn.utils.testing import assert_almost_equal
from sklearn.utils.testing import assert_array_almost_equal
from sklearn.utils.testing import assert_equal
from sklearn.utils.testing import assert_array_equal
from sklearn.utils.testing import assert_greater
from sklearn.utils.testing import assert_raises
from sklearn.utils.testing import assert_raise_message
from sklearn.utils.testing import assert_warns_message
from sklearn.utils.testing import ignore_warnings
from sklearn.utils import check_random_state

from sklearn import datasets
from sklearn.metrics import mean_squared_error
from sklearn.metrics.scorer import SCORERS
from sklearn.metrics import make_scorer

from sklearn.linear_model.base import LinearRegression
from sklearn.linear_model.ridge import ridge_regression
from sklearn.linear_model.ridge import Ridge
from sklearn.linear_model.ridge import _RidgeGCV
from sklearn.linear_model.ridge import RidgeCV
from sklearn.linear_model.ridge import RidgeClassifier
from sklearn.linear_model.ridge import RidgeClassifierCV
from sklearn.linear_model.ridge import _solve_cholesky
from sklearn.linear_model.ridge import _solve_cholesky_kernel
from sklearn.linear_model.ridge import ridge_path
from sklearn.linear_model.ridge import _precomp_kernel_ridge_path_eigen
from sklearn.linear_model.ridge import _kernel_ridge_path_eigen
from sklearn.linear_model.ridge import _ridge_gcv_path_svd

from sklearn.cross_validation import KFold
from sklearn.cross_validation import LeaveOneOut

diabetes = datasets.load_diabetes()
X_diabetes, y_diabetes = diabetes.data, diabetes.target
ind = np.arange(X_diabetes.shape[0])
rng = np.random.RandomState(0)
rng.shuffle(ind)
ind = ind[:200]
X_diabetes, y_diabetes = X_diabetes[ind], y_diabetes[ind]

iris = datasets.load_iris()

X_iris = sp.csr_matrix(iris.data)
y_iris = iris.target

DENSE_FILTER = lambda X: X
SPARSE_FILTER = lambda X: sp.csr_matrix(X)


def test_ridge():
    """Ridge regression convergence test using score

    TODO: for this test to be robust, we should use a dataset instead
    of np.random.
    """
    rng = np.random.RandomState(0)
    alpha = 1.0

    for solver in ("svd", "sparse_cg", "cholesky", "lsqr", "eigen"):
        # With more samples than features
        n_samples, n_features = 6, 5
        y = rng.randn(n_samples)
        X = rng.randn(n_samples, n_features)

        ridge = Ridge(alpha=alpha, solver=solver)
        ridge.fit(X, y)
        assert_equal(ridge.coef_.shape, (X.shape[1], ))
        assert_greater(ridge.score(X, y), 0.47)

        if solver == "cholesky":
            # Currently the only solver to support sample_weight.
            ridge.fit(X, y, sample_weight=np.ones(n_samples))
            assert_greater(ridge.score(X, y), 0.47)

        # With more features than samples
        n_samples, n_features = 5, 10
        y = rng.randn(n_samples)
        X = rng.randn(n_samples, n_features)
        ridge = Ridge(alpha=alpha, solver=solver)
        ridge.fit(X, y)
        assert_greater(ridge.score(X, y), .9)

        if solver == "cholesky":
            # Currently the only solver to support sample_weight.
            ridge.fit(X, y, sample_weight=np.ones(n_samples))
            assert_greater(ridge.score(X, y), 0.9)


def test_primal_dual_relationship():
    y = y_diabetes.reshape(-1, 1)
    coef = _solve_cholesky(X_diabetes, y, alpha=[1e-2])
    K = np.dot(X_diabetes, X_diabetes.T)
    dual_coef = _solve_cholesky_kernel(K, y, alpha=[1e-2])
    coef2 = np.dot(X_diabetes.T, dual_coef).T
    assert_array_almost_equal(coef, coef2)


def test_ridge_singular():
    # test on a singular matrix
    rng = np.random.RandomState(0)
    n_samples, n_features = 6, 6
    y = rng.randn(n_samples // 2)
    y = np.concatenate((y, y))
    X = rng.randn(n_samples // 2, n_features)
    X = np.concatenate((X, X), axis=0)

    ridge = Ridge(alpha=0)
    ridge.fit(X, y)
    assert_greater(ridge.score(X, y), 0.9)


def test_ridge_sample_weights():
    rng = np.random.RandomState(0)

    for solver in ("cholesky", ):
        for n_samples, n_features in ((6, 5), (5, 10)):
            for alpha in (1.0, 1e-2):
                y = rng.randn(n_samples)
                X = rng.randn(n_samples, n_features)
                sample_weight = 1 + rng.rand(n_samples)

                coefs = ridge_regression(X, y,
                                         alpha=alpha,
                                         sample_weight=sample_weight,
                                         solver=solver)
                # Sample weight can be implemented via a simple rescaling
                # for the square loss.
                coefs2 = ridge_regression(
                    X * np.sqrt(sample_weight)[:, np.newaxis],
                    y * np.sqrt(sample_weight),
                    alpha=alpha, solver=solver)
                assert_array_almost_equal(coefs, coefs2)

                # Test for fit_intercept = True
                est = Ridge(alpha=alpha, solver=solver)
                est.fit(X, y, sample_weight=sample_weight)

                # Check using Newton's Method
                # Quadratic function should be solved in a single step.
                # Initialize
                sample_weight = np.sqrt(sample_weight)
                X_weighted = sample_weight[:, np.newaxis] * (
                    np.column_stack((np.ones(n_samples), X)))
                y_weighted = y * sample_weight

                # Gradient is (X*coef-y)*X + alpha*coef_[1:]
                # Remove coef since it is initialized to zero.
                grad = -np.dot(y_weighted, X_weighted)

                # Hessian is (X.T*X) + alpha*I except that the first
                # diagonal element should be zero, since there is no
                # penalization of intercept.
                diag = alpha * np.ones(n_features + 1)
                diag[0] = 0.
                hess = np.dot(X_weighted.T, X_weighted)
                hess.flat[::n_features + 2] += diag
                coef_ = - np.dot(linalg.inv(hess), grad)
                assert_almost_equal(coef_[0], est.intercept_)
                assert_array_almost_equal(coef_[1:], est.coef_)


def test_ridge_shapes():
    """Test shape of coef_ and intercept_
    """
    rng = np.random.RandomState(0)
    n_samples, n_features = 5, 10
    X = rng.randn(n_samples, n_features)
    y = rng.randn(n_samples)
    Y1 = y[:, np.newaxis]
    Y = np.c_[y, 1 + y]

    ridge = Ridge()

    ridge.fit(X, y)
    assert_equal(ridge.coef_.shape, (n_features,))
    assert_equal(ridge.intercept_.shape, ())

    ridge.fit(X, Y1)
    assert_equal(ridge.coef_.shape, (1, n_features))
    assert_equal(ridge.intercept_.shape, (1, ))

    ridge.fit(X, Y)
    assert_equal(ridge.coef_.shape, (2, n_features))
    assert_equal(ridge.intercept_.shape, (2, ))


def test_ridge_intercept():
    """Test intercept with multiple targets GH issue #708
    """
    rng = np.random.RandomState(0)
    n_samples, n_features = 5, 10
    X = rng.randn(n_samples, n_features)
    y = rng.randn(n_samples)
    Y = np.c_[y, 1. + y]

    ridge = Ridge()

    ridge.fit(X, y)
    intercept = ridge.intercept_

    ridge.fit(X, Y)
    assert_almost_equal(ridge.intercept_[0], intercept)
    assert_almost_equal(ridge.intercept_[1], intercept + 1.)


def test_toy_ridge_object():
    """Test BayesianRegression ridge classifier

    TODO: test also n_samples > n_features
    """
    X = np.array([[1], [2]])
    Y = np.array([1, 2])
    clf = Ridge(alpha=0.0)
    clf.fit(X, Y)
    X_test = [[1], [2], [3], [4]]
    assert_almost_equal(clf.predict(X_test), [1., 2, 3, 4])

    assert_equal(len(clf.coef_.shape), 1)
    assert_equal(type(clf.intercept_), np.float64)

    Y = np.vstack((Y, Y)).T

    clf.fit(X, Y)
    X_test = [[1], [2], [3], [4]]

    assert_equal(len(clf.coef_.shape), 2)
    assert_equal(type(clf.intercept_), np.ndarray)


def test_ridge_vs_lstsq():
    """On alpha=0., Ridge and OLS yield the same solution."""

    rng = np.random.RandomState(0)
    # we need more samples than features
    n_samples, n_features = 5, 4
    y = rng.randn(n_samples)
    X = rng.randn(n_samples, n_features)

    ridge = Ridge(alpha=0., fit_intercept=False)
    ols = LinearRegression(fit_intercept=False)

    ridge.fit(X, y)
    ols.fit(X, y)
    assert_almost_equal(ridge.coef_, ols.coef_)

    ridge.fit(X, y)
    ols.fit(X, y)
    assert_almost_equal(ridge.coef_, ols.coef_)


def test_ridge_individual_penalties():
    """Tests the ridge object using individual penalties"""

    rng = np.random.RandomState(42)

    n_samples, n_features, n_targets = 20, 10, 5
    X = rng.randn(n_samples, n_features)
    y = rng.randn(n_samples, n_targets)

    penalties = np.arange(n_targets)

    coef_cholesky = np.array([
        Ridge(alpha=alpha, solver="cholesky").fit(X, target).coef_
        for alpha, target in zip(penalties, y.T)])

    coefs_indiv_pen = [
        Ridge(alpha=penalties, solver=solver, tol=1e-6).fit(X, y).coef_
        for solver in ['svd', 'sparse_cg', 'lsqr', 'cholesky', 'eigen']]
    for coef_indiv_pen in coefs_indiv_pen:
        assert_array_almost_equal(coef_cholesky, coef_indiv_pen)

    # Test error is raised when number of targets and penalties do not match.
    ridge = Ridge(alpha=penalties[:3])
    assert_raises(ValueError, ridge.fit, X, y)


def _test_ridge_loo(filter_):
    # test that can work with both dense or sparse matrices
    n_samples = X_diabetes.shape[0]

    ret = []

    ridge_gcv = _RidgeGCV(fit_intercept=False)
    ridge = Ridge(alpha=1.0, fit_intercept=False)

    alphas = np.array([1.])
    errors, _ = _kernel_ridge_path_eigen(X_diabetes,
                                         y_diabetes[:, np.newaxis],
                                         alphas,
                                         mode='looe')
    errors = errors.ravel() ** 2
    values, _ = _kernel_ridge_path_eigen(X_diabetes,
                                         y_diabetes[:, np.newaxis],
                                         alphas,
                                         mode='loov')
    values = values.ravel()

    # brute-force leave-one-out: remove one example at a time
    errors2 = []
    values2 = []
    for i in range(n_samples):
        sel = np.arange(n_samples) != i
        X_new = X_diabetes[sel]
        y_new = y_diabetes[sel]
        ridge.fit(X_new, y_new)
        value = ridge.predict([X_diabetes[i]])[0]
        error = (y_diabetes[i] - value) ** 2
        errors2.append(error)
        values2.append(value)

    # check that efficient and brute-force LOO give same results
    assert_almost_equal(errors, errors2)
    assert_almost_equal(values, values2)

    # generalized cross-validation (efficient leave-one-out,
    # SVD variation)
    errors3_, c = _ridge_gcv_path_svd(X_diabetes,
                                      np.atleast_2d(y_diabetes.T).T,
                                      np.atleast_2d(ridge.alpha),
                                      mode='looe')
    errors3 = (errors3_ ** 2).ravel()
    values3, c = _ridge_gcv_path_svd(X_diabetes,
                                      np.atleast_2d(y_diabetes.T).T,
                                      np.atleast_2d(ridge.alpha),
                                      mode='loov')
    values3 = values3.ravel()

    # check that efficient and SVD efficient LOO give same results
    assert_almost_equal(errors, errors3)
    assert_almost_equal(values, values3)

    # check best alpha
    ridge_gcv.fit(filter_(X_diabetes), y_diabetes)
    alpha_ = ridge_gcv.alpha_
    ret.append(alpha_)

    # check that we get same best alpha with custom loss_func
    f = ignore_warnings
    scoring = make_scorer(mean_squared_error, greater_is_better=False)
    ridge_gcv2 = RidgeCV(fit_intercept=False, scoring=scoring)
    f(ridge_gcv2.fit)(filter_(X_diabetes), y_diabetes)
    assert_equal(ridge_gcv2.alpha_, alpha_)

    # check that we get same best alpha with custom score_func
    func = lambda x, y: -mean_squared_error(x, y)
    scoring = make_scorer(func)
    ridge_gcv3 = RidgeCV(fit_intercept=False, scoring=scoring)
    f(ridge_gcv3.fit)(filter_(X_diabetes), y_diabetes)
    assert_equal(ridge_gcv3.alpha_, alpha_)

    # check that we get same best alpha with a scorer
    scorer = SCORERS['mean_squared_error']
    ridge_gcv4 = RidgeCV(fit_intercept=False, scoring=scorer)
    ridge_gcv4.fit(filter_(X_diabetes), y_diabetes)
    assert_equal(ridge_gcv4.alpha_, alpha_)

    # check that we get same best alpha with sample weights
    ridge_gcv.fit(filter_(X_diabetes), y_diabetes,
                  sample_weight=np.ones(n_samples))
    assert_equal(ridge_gcv.alpha_, alpha_)

    # simulate several responses
    Y = np.vstack((y_diabetes, y_diabetes)).T

    ridge_gcv.fit(filter_(X_diabetes), Y)
    Y_pred = ridge_gcv.predict(filter_(X_diabetes))
    ridge_gcv.fit(filter_(X_diabetes), y_diabetes)
    y_pred = ridge_gcv.predict(filter_(X_diabetes))

    assert_array_almost_equal(np.vstack((y_pred, y_pred)).T,
                              Y_pred, decimal=5)

    return ret


def _test_ridge_cv(filter_):
    n_samples = X_diabetes.shape[0]

    ridge_cv = RidgeCV()
    ridge_cv.fit(filter_(X_diabetes), y_diabetes)
    ridge_cv.predict(filter_(X_diabetes))

    assert_equal(len(ridge_cv.coef_.shape), 1)
    assert_equal(type(ridge_cv.intercept_), np.float64)

    cv = KFold(n_samples, 5)
    ridge_cv.set_params(cv=cv)
    ridge_cv.fit(filter_(X_diabetes), y_diabetes)
    ridge_cv.predict(filter_(X_diabetes))

    assert_equal(len(ridge_cv.coef_.shape), 1)
    assert_equal(type(ridge_cv.intercept_), np.float64)


def _test_ridge_diabetes(filter_):
    ridge = Ridge(fit_intercept=False)
    ridge.fit(filter_(X_diabetes), y_diabetes)
    return np.round(ridge.score(filter_(X_diabetes), y_diabetes), 5)


def _test_multi_ridge_diabetes(filter_):
    # simulate several responses
    Y = np.vstack((y_diabetes, y_diabetes)).T
    n_features = X_diabetes.shape[1]

    ridge = Ridge(fit_intercept=False)
    ridge.fit(filter_(X_diabetes), Y)
    assert_equal(ridge.coef_.shape, (2, n_features))
    Y_pred = ridge.predict(filter_(X_diabetes))
    ridge.fit(filter_(X_diabetes), y_diabetes)
    y_pred = ridge.predict(filter_(X_diabetes))
    assert_array_almost_equal(np.vstack((y_pred, y_pred)).T,
                              Y_pred, decimal=3)


def _test_ridge_classifiers(filter_):
    n_classes = np.unique(y_iris).shape[0]
    n_features = X_iris.shape[1]
    for clf in (RidgeClassifier(), RidgeClassifierCV()):
        clf.fit(filter_(X_iris), y_iris)
        assert_equal(clf.coef_.shape, (n_classes, n_features))
        y_pred = clf.predict(filter_(X_iris))
        assert_greater(np.mean(y_iris == y_pred), .79)

    n_samples = X_iris.shape[0]
    cv = KFold(n_samples, 5)
    clf = RidgeClassifierCV(cv=cv)
    clf.fit(filter_(X_iris), y_iris)
    y_pred = clf.predict(filter_(X_iris))
    assert_true(np.mean(y_iris == y_pred) >= 0.8)


def _test_tolerance(filter_):
    ridge = Ridge(tol=1e-5)
    ridge.fit(filter_(X_diabetes), y_diabetes)
    score = ridge.score(filter_(X_diabetes), y_diabetes)

    ridge2 = Ridge(tol=1e-3)
    ridge2.fit(filter_(X_diabetes), y_diabetes)
    score2 = ridge2.score(filter_(X_diabetes), y_diabetes)

    assert_true(score >= score2)


def test_dense_sparse():
    for test_func in (_test_ridge_loo,
                      _test_ridge_cv,
                      _test_ridge_diabetes,
                      _test_multi_ridge_diabetes,
                      _test_ridge_classifiers,
                      _test_tolerance):
        # test dense matrix
        ret_dense = test_func(DENSE_FILTER)
        # test sparse matrix
        ret_sparse = test_func(SPARSE_FILTER)
        # test that the outputs are the same
        if ret_dense is not None and ret_sparse is not None:
            assert_array_almost_equal(ret_dense, ret_sparse, decimal=3)


def test_ridge_cv_sparse_svd():
    X = sp.csr_matrix(X_diabetes)
    ridge = RidgeCV(gcv_mode="svd")
    assert_raises(TypeError, ridge.fit, X)


def test_class_weights():
    """
    Test class weights.
    """
    X = np.array([[-1.0, -1.0], [-1.0, 0], [-.8, -1.0],
                  [1.0, 1.0], [1.0, 0.0]])
    y = [1, 1, 1, -1, -1]

    clf = RidgeClassifier(class_weight=None)
    clf.fit(X, y)
    assert_array_equal(clf.predict([[0.2, -1.0]]), np.array([1]))

    # we give a small weights to class 1
    clf = RidgeClassifier(class_weight={1: 0.001})
    clf.fit(X, y)

    # now the hyperplane should rotate clock-wise and
    # the prediction on this point should shift
    assert_array_equal(clf.predict([[0.2, -1.0]]), np.array([-1]))

    # check if class_weight = 'auto' can handle negative labels.
    clf = RidgeClassifier(class_weight='auto')
    clf.fit(X, y)
    assert_array_equal(clf.predict([[0.2, -1.0]]), np.array([1]))

    # class_weight = 'auto', and class_weight = None should return
    # same values when y has equal number of all labels
    X = np.array([[-1.0, -1.0], [-1.0, 0], [-.8, -1.0], [1.0, 1.0]])
    y = [1, 1, -1, -1]
    clf = RidgeClassifier(class_weight=None)
    clf.fit(X, y)
    clfa = RidgeClassifier(class_weight='auto')
    clfa.fit(X, y)
    assert_equal(len(clfa.classes_), 2)
    assert_array_almost_equal(clf.coef_, clfa.coef_)
    assert_array_almost_equal(clf.intercept_, clfa.intercept_)


def test_class_weights_cv():
    """
    Test class weights for cross validated ridge classifier.
    """
    X = np.array([[-1.0, -1.0], [-1.0, 0], [-.8, -1.0],
                  [1.0, 1.0], [1.0, 0.0]])
    y = [1, 1, 1, -1, -1]

    clf = RidgeClassifierCV(class_weight=None, alphas=[.01, .1, 1])
    clf.fit(X, y)

    # we give a small weights to class 1
    clf = RidgeClassifierCV(class_weight={1: 0.001}, alphas=[.01, .1, 1, 10])
    clf.fit(X, y)

    assert_array_equal(clf.predict([[-.2, 2]]), np.array([-1]))


def test_ridgecv_store_cv_values():
    """
    Test _RidgeCV's store_cv_values attribute.
    """
    rng = rng = np.random.RandomState(42)

    n_samples = 8
    n_features = 5
    x = rng.randn(n_samples, n_features)
    alphas = [1e-1, 1e0, 1e1]
    n_alphas = len(alphas)

    r = RidgeCV(alphas=alphas, store_cv_values=True)

    # with len(y.shape) == 1
    y = rng.randn(n_samples)
    r.fit(x, y)
    assert_equal(r.cv_values_.shape, (n_samples, n_alphas))

    # with len(y.shape) == 2
    n_responses = 3
    y = rng.randn(n_samples, n_responses)
    r.fit(x, y)
    assert_equal(r.cv_values_.shape, (n_samples, n_responses, n_alphas))


def test_ridge_sample_weights_in_feature_space():
    """Check that Cholesky solver in feature space applies sample_weights
    correctly.
    """

    rng = np.random.RandomState(42)

    n_samples_list = [5, 6, 7] * 2
    n_features_list = [7, 6, 5] * 2
    n_targets_list = [1, 1, 1, 2, 2, 2]
    noise = 1.
    alpha = 2.
    alpha = np.atleast_1d(alpha)

    for n_samples, n_features, n_targets in zip(n_samples_list,
                                                n_features_list,
                                                n_targets_list):
        X = rng.randn(n_samples, n_features)
        beta = rng.randn(n_features, n_targets)
        Y = X.dot(beta)
        Y_noisy = Y + rng.randn(*Y.shape) * np.sqrt((Y ** 2).sum(0)) * noise

        K = X.dot(X.T)
        sample_weights = 1. + (rng.randn(n_samples) ** 2) * 10

        coef_sample_space = _solve_cholesky_kernel(K, Y_noisy, alpha,
                                         sample_weight=sample_weights)

        coef_feature_space = _solve_cholesky(X, Y_noisy, alpha,
                                         sample_weight=sample_weights)

        assert_array_almost_equal(X.T.dot(coef_sample_space),
                                  coef_feature_space.T)


def test_raises_value_error_if_sample_weights_greater_than_1d():
    """Sample weights must be either scalar or 1D"""

    n_sampless = [2, 3]
    n_featuress = [3, 2]

    rng = np.random.RandomState(42)


    for n_samples, n_features in zip(n_sampless, n_featuress):
        X = rng.randn(n_samples, n_features)
        y = rng.randn(n_samples)
        sample_weights_OK = rng.randn(n_samples) ** 2 + 1
        sample_weights_OK_1 = 1.
        sample_weights_OK_2 = 2.
        sample_weights_not_OK = sample_weights_OK[:, np.newaxis]
        sample_weights_not_OK_2 = sample_weights_OK[np.newaxis, :]

        ridge = Ridge(alpha=1)

        # make sure the "OK" sample weights actually work
        ridge.fit(X, y, sample_weights_OK)
        ridge.fit(X, y, sample_weights_OK_1)
        ridge.fit(X, y, sample_weights_OK_2)

        def fit_ridge_not_ok():
            ridge.fit(X, y, sample_weights_not_OK)

        def fit_ridge_not_ok_2():
            ridge.fit(X, y, sample_weights_not_OK_2)

        assert_raise_message(ValueError,
                              "Sample weights must be 1D array or scalar",
                              fit_ridge_not_ok)

        assert_raise_message(ValueError,
                              "Sample weights must be 1D array or scalar",
                              fit_ridge_not_ok_2)


def test_sparse_design_with_sample_weights():
    """Sample weights must work with sparse matrices"""

    n_sampless = [2, 3]
    n_featuress = [3, 2]

    rng = np.random.RandomState(42)

    sparse_matrix_converters = [sp.coo_matrix,
                                sp.csr_matrix,
                                sp.csc_matrix,
                                sp.lil_matrix,
                                sp.dok_matrix
                                ]

    sparse_ridge = Ridge(alpha=1., fit_intercept=False)
    dense_ridge = Ridge(alpha=1., fit_intercept=False)

    for n_samples, n_features in zip(n_sampless, n_featuress):
        X = rng.randn(n_samples, n_features)
        y = rng.randn(n_samples)
        sample_weights = rng.randn(n_samples) ** 2 + 1
        for sparse_converter in sparse_matrix_converters:
            X_sparse = sparse_converter(X)
            sparse_ridge.fit(X_sparse, y, sample_weight=sample_weights)
            dense_ridge.fit(X, y, sample_weight=sample_weights)

            assert_array_almost_equal(sparse_ridge.coef_, dense_ridge.coef_,
                                      decimal=6)


def test_deprecation_warning_dense_cholesky():
    """Tests if DeprecationWarning is raised at instantiation of estimators
    and when ridge_regression is called"""

    warning_class = DeprecationWarning
    warning_message = ("The name 'dense_cholesky' is deprecated."
                       " Using 'cholesky' instead")

    X = np.ones([2, 3])
    y = np.ones(2)
    func1 = lambda: Ridge(solver='dense_cholesky').fit(X, y)
    func2 = lambda: RidgeClassifier(solver='dense_cholesky').fit(X, y)
    X = np.ones([3, 2])
    y = np.zeros(3)
    func3 = lambda: ridge_regression(X, y, alpha=1, solver='dense_cholesky')

    for func in [func1, func2, func3]:
        assert_warns_message(warning_class, warning_message, func)


def test_raises_value_error_if_solver_not_supported():
    """Tests whether a ValueError is raised if a non-identified solver
    is passed to ridge_regression"""

    wrong_solver = "This is not a solver (MagritteSolveCV QuantumBitcoin)"

    exception = ValueError
    message = "Solver %s not understood" % wrong_solver

    def func():
        X = np.eye(3)
        y = np.ones(3)
        ridge_regression(X, y, alpha=1., solver=wrong_solver)

    assert_raise_message(exception, message, func)


def test_ridge_cv_path_solvers():

    alphas = np.logspace(-3, 3, 10)
    cv = KFold(len(y_diabetes), 2)

    for solver in ['eigen', 'svd']:
        ridge_cv = RidgeCV(alphas=alphas, solver=None, cv=cv,
                       fit_intercept=False)
        ridge_cv_path = RidgeCV(alphas=alphas, solver=solver, cv=cv,
                                fit_intercept=False)

        ridge_cv.fit(X_diabetes, y_diabetes)
        ridge_cv_path.fit(X_diabetes, y_diabetes[:, np.newaxis])
        assert_array_almost_equal(ridge_cv.coef_, ridge_cv_path.coef_[0])

        ridge_cv.fit_intercept = True
        ridge_cv_path.fit_intercept = True

        ridge_cv.fit(X_diabetes, y_diabetes)
        ridge_cv_path.fit(X_diabetes, y_diabetes[:, np.newaxis])

        assert_array_almost_equal(ridge_cv.coef_, ridge_cv_path.coef_[0])


def make_noisy_forward_data(
    n_samples=100,
    n_features=200,
    n_targets=10,
    train_frac=.8,
    noise_levels=None,
    random_state=42):
    """Creates a simple, dense, noisy forward linear model with multiple
    output."""
    rng = check_random_state(random_state)
    n_train = int(train_frac * n_samples)
    train = slice(None, n_train)
    test = slice(n_train, None)
    X = rng.randn(n_samples, n_features)
    W = rng.randn(n_features, n_targets)
    Y_clean = X.dot(W)
    if noise_levels is None:
        noise_levels = rng.randn(n_targets) ** 2
    noise_levels = np.atleast_1d(noise_levels) * np.ones(n_targets)
    noise = rng.randn(*Y_clean.shape) * noise_levels * Y_clean.std(0)
    Y = Y_clean + noise
    return X, Y, W, train, test


def test_ridge_path():

    n_sampless, n_featuress, n_targetss = [100, 200], [200, 100], [10, 10]
    for n_samples, n_features, n_targets in zip(n_sampless,
                                                n_featuress,
                                                n_targetss):
        X, Y, W, train, test = make_noisy_forward_data(
            n_samples, n_features, n_targets)
        alphas = np.logspace(-3, 3, 9)[:, np.newaxis]

        predictions_cholesky = np.array([Ridge(
                    solver='cholesky', alpha=alpha, fit_intercept=False).fit(
                    X[train], Y[train]).predict(X[test])
                                         for alpha in alphas])

        predictions_eigen_path = ridge_path(X[train], Y[train], alphas,
                                            X[test], solver='eigen')

        predictions_svd_path = ridge_path(X[train], Y[train], alphas,
                                          X[test], solver='svd')

        assert_array_almost_equal(predictions_cholesky,
                                  predictions_eigen_path)
        assert_array_almost_equal(predictions_cholesky,
                                  predictions_svd_path)


def test__precomp_kernel_ridge_path_eigen_loov():

    n_samples, n_features, n_targets = 20, 40, 10
    X, Y, W, _, _ = make_noisy_forward_data(
            n_samples, n_features, n_targets)

    alphas = np.logspace(-3, 3, 9)
    gramX = X.dot(X.T)
    cv = LeaveOneOut(n_samples)

    loov = _precomp_kernel_ridge_path_eigen(gramX, Y, alphas[:, np.newaxis],
                                            mode='loov')[0]

    loov_normal = [
        _precomp_kernel_ridge_path_eigen(
            gramX[train[:, np.newaxis], train], Y[train], 
            alphas[:, np.newaxis],
            gramX[test[:, np.newaxis], train])
        for train, test in cv]

    assert_array_almost_equal(
        loov, np.array(loov_normal).squeeze().transpose(1, 0, 2))


def test_ridge_svd_tall():
    """Make sure the SVD approach also works for tall matrices"""

    n_samples, n_features, n_targets = 20, 5, 2
    X, Y, W, _, _ = make_noisy_forward_data(
        n_samples, n_features, n_targets)

    coef_svd = Ridge(alpha=1, solver='svd').fit(X, Y).coef_
    coef_cholesky = Ridge(alpha=1, solver='cholesky').fit(X, Y).coef_

    assert_array_almost_equal(coef_svd, coef_cholesky)


def test_ridge_gcv_looe_on_dense_data():

    n_samples, n_features, n_targets = 50, 10, 2
    X, Y, W, _, _ = make_noisy_forward_data(n_samples, n_features, n_targets)

    alphas = np.logspace(-3, 3, 9)

    sample_weights = np.ones(len(X))

    # adding sample weights will go to the old branch, without it will go
    # to the new branch

    old_coef = _RidgeGCV(alphas=alphas, gcv_mode="eigen").fit(
        X, Y, sample_weight=sample_weights).coef_

    new_coef = _RidgeGCV(alphas=alphas, gcv_mode="eigen").fit(X, Y).coef_

    assert_array_almost_equal(old_coef, new_coef)


def test_ridge_gcv_with_sample_weights():

    n_samples, n_features, n_targets = 20, 5, 7
    X, Y, W, _, _ = make_noisy_forward_data(n_samples, n_features, n_targets)
    alphas = np.logspace(-3, 3, 9)

    rng = np.random.RandomState(42)
    sample_weights = rng.randn(n_samples) ** 2

    cv = LeaveOneOut(n_samples)
    cv_predictions = np.array([[
        Ridge(solver='cholesky', alpha=alpha, fit_intercept=False).fit(
            X[train], Y[train], sample_weight=sample_weights[train]
            ).predict(X[test])
        for train, test in cv] for alpha in alphas]).squeeze()

    cv_errors = Y[np.newaxis] - cv_predictions.reshape(
        len(alphas), n_samples, n_targets)

    ridge_gcv = _RidgeGCV(alphas=alphas, store_cv_values=True,
                          gcv_mode='eigen', fit_intercept=False)
    # emulate the sample weight stuff from _RidgeGCV
    ridge_gcv.fit(X, Y, sample_weight=sample_weights)
    loo_predictions = ridge_gcv.cv_values_


    assert_array_almost_equal(cv_errors ** 2,
                              loo_predictions.transpose(2, 0, 1))


def test_kernel_ridge_path_with_sample_weights():
    n_samples, n_features, n_targets = 20, 5, 7
    X, Y, W, _, _ = make_noisy_forward_data(n_samples, n_features, n_targets)
    alphas = np.logspace(-3, 3, 9)[:, np.newaxis] * \
        np.arange(1, n_targets + 1)

    rng = np.random.RandomState(42)
    sample_weights = rng.randn(n_samples) ** 2

    cv = LeaveOneOut(n_samples)
    cv_predictions = np.array([[
        Ridge(solver='cholesky', alpha=alpha, fit_intercept=False).fit(
            X[train], Y[train], sample_weight=sample_weights[train]
            ).predict(X[test])
        for train, test in cv] for alpha in alphas]).squeeze()

    cv_errors = Y[np.newaxis] - cv_predictions.reshape(
        len(alphas), n_samples, n_targets)

    ridge_path_gcv_errors = _kernel_ridge_path_eigen(
        X, Y, alphas, sample_weight=sample_weights, mode='looe')[0]

    assert_array_almost_equal(cv_errors, ridge_path_gcv_errors)


def test__ridge_gcv_path_svd_against_eigen():
    n_samples, n_features, n_targets = 20, 5, 7
    X, Y, W, _, _ = make_noisy_forward_data(n_samples, n_features, n_targets)
    alphas = np.logspace(-3, 3, 9)[:, np.newaxis] * \
        np.arange(1, n_targets + 1)

    eigen_solved, c = _kernel_ridge_path_eigen(X, Y, alphas, mode='looe')
    svd_solved, c2 = _ridge_gcv_path_svd(X, Y, alphas, mode='looe')


    assert_array_almost_equal(eigen_solved, svd_solved)


def test__ridge_gcv_path_svd_with_sample_weights_against_eigen():
    n_samples, n_features, n_targets = 20, 5, 7
    X, Y, W, _, _ = make_noisy_forward_data(n_samples, n_features, n_targets)
    alphas = np.logspace(-3, 3, 9)[:, np.newaxis] * \
        np.arange(1, n_targets + 1)

    rng = np.random.RandomState(42)
    sample_weights = rng.randn(n_samples) ** 2

    looe_eigen = _kernel_ridge_path_eigen(
        X, Y, alphas, sample_weight=sample_weights, mode='looe')[0]
    looe_svd = _ridge_gcv_path_svd(X, Y, alphas,
                                   sample_weight=sample_weights,
                                   mode='looe')[0]

    assert_array_almost_equal(looe_eigen, looe_svd)


def test_best_alpha_scales_with_target():
    """Best penalty must be selected per target and yield corresponding
    coefficients"""

    n_samples, n_features, n_targets = 20, 10, 10
    noise_levels = np.arange(10) / 5.
    X, Y, W, _, _ = make_noisy_forward_data(
        n_samples, n_features, n_targets, noise_levels=noise_levels)

    alphas = np.logspace(-3, 3, 35)

    cv_mses = []
    coefs = []
    best_alphas = []
    for y in Y.T:
        ridge_gcv = _RidgeGCV(alphas=alphas, store_cv_values=True,
                          fit_intercept=False)

        ridge_gcv.fit(X, y)
        cv_mses.append(ridge_gcv.cv_values_.mean(0))
        coefs.append(ridge_gcv.coef_)
        best_alphas.append(ridge_gcv.alpha_)

    ridge_gcv.fit(X, Y)
    assert_array_almost_equal(
        np.array(cv_mses), ridge_gcv.cv_values_.mean(0))
    assert_array_almost_equal(np.array(coefs),
                              ridge_gcv.coef_)



