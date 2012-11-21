import numpy as np
import scipy.sparse as sp

from sklearn.utils.testing import assert_true
from sklearn.utils.testing import assert_almost_equal
from sklearn.utils.testing import assert_array_almost_equal
from sklearn.utils.testing import assert_equal
from sklearn.utils.testing import assert_array_equal
from sklearn.utils.testing import assert_greater

from sklearn import datasets
from sklearn.metrics import mean_squared_error

from sklearn.linear_model.base import LinearRegression
from sklearn.linear_model.ridge import Ridge
from sklearn.linear_model.ridge import _RidgeGCV
from sklearn.linear_model.ridge import RidgeCV
from sklearn.linear_model.ridge import RidgeClassifier
from sklearn.linear_model.ridge import RidgeClassifierCV
from sklearn.linear_model.ridge import ridge_regression
from sklearn.linear_model.ridge import _RidgeGridCV

from sklearn.cross_validation import KFold
from sklearn.metrics import euclidean_distances
from sklearn.datasets import make_regression

from itertools import product
import numbers

rng = np.random.RandomState(0)
diabetes = datasets.load_diabetes()
X_diabetes, y_diabetes = diabetes.data, diabetes.target
ind = np.arange(X_diabetes.shape[0])
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
    alpha = 1.0

    for solver in ("sparse_cg", "dense_cholesky", "lsqr", "svd", "eigen"):
        # With more samples than features
        n_samples, n_features = 6, 5
        y = rng.randn(n_samples)
        X = rng.randn(n_samples, n_features)

        ridge = Ridge(alpha=alpha, solver=solver)
        ridge.fit(X, y)
        assert_equal(ridge.coef_.shape, (X.shape[1], ))
        assert_greater(ridge.score(X, y), 0.47)

        ridge.fit(X, y, sample_weight=np.ones(n_samples))
        assert_greater(ridge.score(X, y), 0.47)

        # With more features than samples
        n_samples, n_features = 5, 10
        y = rng.randn(n_samples)
        X = rng.randn(n_samples, n_features)
        ridge = Ridge(alpha=alpha, solver=solver)
        ridge.fit(X, y)
        assert_greater(ridge.score(X, y), .9)

        ridge.fit(X, y, sample_weight=np.ones(n_samples))
        assert_greater(ridge.score(X, y), 0.9)


def test_ridge_compare_different_solvers():
    """Compares regression results for all solvers

    This test covers the cases n_samples > n_features
    and n_samples < n_features using different penalties
    on a noisy linear model with exactly one target.
    """

    tolerance = 0.02
    # sparse_cg solver differs significantly from the others. Is this
    # tolerance too permissive? (We are checking for general correctness
    # of solution, nothing more)
    # Somehow even with a fixed random seed in the data, sparse_cg seems
    # to output different results each time. I guess this is due to design
    # (randomness in gradient descent?). Set tolerance to e.g 0.007 and the
    # tests will not pass every time. (Even 0.02 doesn't pass every time)

    alphas = [0., .01, .1, 1., 10.]

    solvers = ["sparse_cg", "dense_cholesky", "lsqr", "svd", "eigen"]

    N, P = 20, 50
    noise_level = .1

    # N = n_samples < n_features = P
    X_1 = rng.randn(N, P)
    beta_1 = rng.randn(P)
    noise_1 = noise_level * rng.randn(N)
    y_1 = np.dot(X_1, beta_1) + noise_1

    # P = n_samples > n_features = N
    X_2 = rng.randn(P, N)
    beta_2 = rng.randn(N)
    noise_2 = noise_level * rng.randn(P)
    y_2 = np.dot(X_2, beta_2) + noise_2

    test_cases = [(X_1, y_1), (X_2, y_2)]

    for i, (X, y) in enumerate(test_cases):
        for alpha in alphas:
            coef_vectors = [ridge_regression(X, y, alpha, solver=solver)
                            for solver in solvers]

            # pairwise comparison using sklearn.metric.euclidean_distances
            coef_vectors = np.array(coef_vectors)
            distances = euclidean_distances(coef_vectors)

            # make sure these distances are close to 0
            assert_greater(tolerance, distances.max())


# TODO: This test doesn't check dense_cholesky, option multiple targets
# with unique penalty.
def test_ridge_multiple_targets_multiple_penalties():
    """Tests multiple target, multiple individual penalties feature

    against standard cholesky solver applied to each combination
    individually"""

    tolerance = 0.01  # 1e-10 works with svd and eigen
    matrix_shapes = [(20, 50), (50, 20)]
    n_targets = 10
    n_penalties_per_target = 4
    noise_level = .01
    relevant_features_proportion = 0.25

    # for each target generate penalties on a logarithmic scale
    # and multiply them by a target-dependent value to individualize
    lowest_penalty_exponent = 2
    alphas = np.logspace(lowest_penalty_exponent,
                         lowest_penalty_exponent + n_penalties_per_target,
                          n_penalties_per_target, False)
    alphas = alphas[:, np.newaxis] * (rng.rand(1, n_targets) * 9 + 1)

    concerned_solvers = ["svd", "eigen", "lsqr", "sparse_cg", "dense_cholesky"]

    def make_test_case(n_samples, n_features):
        n_informative = int(np.floor(n_features *\
                                     relevant_features_proportion))
        return make_regression(n_samples=n_samples,
                               n_features=n_features,
                               n_informative=n_informative,
                               n_targets=n_targets,
                               noise=noise_level,
                               random_state=rng,
                               coef=False)

    test_cases = [make_test_case(n_samples, n_features)
                  for n_samples, n_features in matrix_shapes]

    # Calculate everything individually using standard solver
    standard_solutions = []
    for X, y in test_cases:
        case_solutions = np.empty([alphas.shape[0],
                                   alphas.shape[1], X.shape[1]])
        for i, alpha_line in enumerate(alphas):
            for j, (target, alpha) in enumerate(zip(y.T, alpha_line)):
                case_solutions[i, j, :] = ridge_regression(
                    X, target, alpha, solver="dense_cholesky")
        standard_solutions.append(case_solutions)

    # now do the same with multiple targets/individual penalties
    new_solutions = []
    for X, y in test_cases:
        case_solutions = np.empty([len(concerned_solvers),
                        alphas.shape[0], alphas.shape[1], X.shape[1]])
        for s, solver in enumerate(concerned_solvers):
            case_solutions[s] = ridge_regression(X, y, alphas, solver=solver)
        new_solutions.append(case_solutions)

    # Compare all these solutions, distances on each target individually
    for standard_solution, new_solution in \
                    zip(standard_solutions, new_solutions):
        distances = (
            (standard_solution[np.newaxis, :] - new_solution) ** 2).sum(-1)

        assert_greater(tolerance, distances.max())


def test_ridge_regression_with_varying_sample_weights():
    pass


def test_ridge_shapes():
    """Test shape of coef_ and intercept_
    """
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


def test_ridge_regression_coef_shapes_individual_penalties():
    """Test shape of coefficients output by ridge_regression
       in the presence of multiple targets/
       multiple individual penalities
    """

    n_sampless = [5, 10]
    n_featuress = [10, 5]
    n_informatives = [2, 2]
    noise_levels = [.1, .1]

    n_targetss = [1, 2, 5]

    datasets = [make_regression(n_samples=n_samples,
                                n_features=n_features,
                                n_informative=n_informative,
                                noise=noise_level,
                                random_state=rng,
                                coef=False,
                                n_targets=n_targets
                                )
                for ((n_samples, n_features, n_informative, noise_level),
                     n_targets) in
                     product(
                        zip(n_sampless, n_featuress,
                             n_informatives, noise_levels),
                        n_targetss)]

    solvers = ["sparse_cg", "lsqr", "svd", "eigen", "dense_cholesky"]

    def expected_shape(X, y, alphas):

        if y.ndim == 1:
            if isinstance(alphas, numbers.Number):
                expected_shape = (X.shape[1],)
            elif alphas.shape == (1,):
                expected_shape = (X.shape[1],)
            elif alphas.ndim == 1:
                expected_shape = (alphas.shape[0], X.shape[1])
            else:
                if alphas.shape[-1] == 1:
                    expected_shape = tuple(list(alphas.shape[:-1]) +
                                           [X.shape[1]])
                else:
                    expected_shape = tuple(list(alphas.shape) + [X.shape[1]])
        else:
            if isinstance(alphas, numbers.Number):
                expected_shape = (y.shape[1], X.shape[1])
            elif alphas.shape[-1] == y.shape[1]:
                expected_shape = tuple(list(alphas.shape) + [X.shape[1]])
            elif alphas.shape[-1] == 1:
                expected_shape = tuple(list(alphas.shape[:-1]) +\
                                    [y.shape[1], X.shape[1]])
            else:
                expected_shape = tuple(list(alphas.shape) +\
                                    [y.shape[1], X.shape[1]])

        return expected_shape

    def make_alphas(y):
        alphas = []

        # one number
        alpha_number = 10.0
        alphas.append(alpha_number)

        # array with one element and one dimension
        alpha_1D_1E_array = np.array([11.0])
        alphas.append(alpha_1D_1E_array)

        # array with one element and two dimensions
        alpha_2D_1E_array = np.array([[12.0]])
        alphas.append(alpha_2D_1E_array)

        # if multiple targets
        if y.ndim == 2:
            # 1D array with len == n_targets
            n_targets = y.shape[1]
            alpha = np.arange(13., 13. + n_targets)
            alphas.append(alpha)

            # 2D array with last dim == n_targets
            scale_factors = np.array([[1.], [2.]])
            alpha2 = scale_factors * alpha[np.newaxis, :]
            alphas.append(alpha2)

            # 1D array with len != n_targets
            alpha3 = np.arange(20., 20. + n_targets + 1)
            alphas.append(alpha3)
        elif y.ndim == 1:
            # 1D array with len != 1
            alpha = np.array([25., 26.])
            alphas.append(alpha)

        # 2D array with last dim == 1
        alpha = np.array([[27.], [28.]])
        alphas.append(alpha)

        # 3D array
        alpha = np.arange(8).reshape(2, 2, 2)
        alphas.append(alpha)

        # 3D array with last dim == 1
        alpha = np.arange(4).reshape(2, 2, 1)
        alphas.append(alpha)

        return alphas

    def verify_shape(X, y, alpha, solver):
        expected = expected_shape(X, y, alpha)
        actual = ridge_regression(X, y, alpha, solver=solver).shape

        assert expected == actual

    for solver in solvers:
        for X, y in datasets:
            alphas = make_alphas(y)
            for alpha in alphas:
                # if one target, check also degenerate 2D y
                if y.ndim == 1:
                    verify_shape(X, y[:, np.newaxis], alpha, solver)
                verify_shape(X, y, alpha, solver)


def test_ridge_intercept():
    """Test intercept with multiple targets GH issue #708
    """
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


def _test_ridge_loo(filter_):
    # test that can work with both dense or sparse matrices
    n_samples = X_diabetes.shape[0]

    ret = []

    ridge_gcv = _RidgeGCV(fit_intercept=False)
    ridge = Ridge(alpha=1.0, fit_intercept=False)

    # generalized cross-validation (efficient leave-one-out)
    decomp = ridge_gcv._pre_compute(X_diabetes, y_diabetes)
    errors, c = ridge_gcv._errors(1.0, y_diabetes, *decomp)
    values, c = ridge_gcv._values(1.0, y_diabetes, *decomp)

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
    decomp = ridge_gcv._pre_compute_svd(X_diabetes, y_diabetes)
    errors3, c = ridge_gcv._errors_svd(ridge.alpha, y_diabetes, *decomp)
    values3, c = ridge_gcv._values_svd(ridge.alpha, y_diabetes, *decomp)

    # check that efficient and SVD efficient LOO give same results
    assert_almost_equal(errors, errors3)
    assert_almost_equal(values, values3)

    # check best alpha
    ridge_gcv.fit(filter_(X_diabetes), y_diabetes)
    alpha_ = ridge_gcv.alpha_
    ret.append(alpha_)

    # check that we get same best alpha with custom loss_func
    ridge_gcv2 = RidgeCV(fit_intercept=False, loss_func=mean_squared_error)
    ridge_gcv2.fit(filter_(X_diabetes), y_diabetes)
    assert_equal(ridge_gcv2.alpha_, alpha_)

    # check that we get same best alpha with custom score_func
    func = lambda x, y: -mean_squared_error(x, y)
    ridge_gcv3 = RidgeCV(fit_intercept=False, score_func=func)
    ridge_gcv3.fit(filter_(X_diabetes), y_diabetes)
    assert_equal(ridge_gcv3.alpha_, alpha_)

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
        if ret_dense != None and ret_sparse != None:
            assert_array_almost_equal(ret_dense, ret_sparse, decimal=3)


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


def test_ridge_grid_cv_object():

    ridge_cv = _RidgeGridCV(n_grid_refinements=10)

    n_samples, n_features = 50, 100

    X = rng.randn(n_samples, n_features)

    beta = np.zeros(n_features)
    beta[:2] = rng.randn(2) + 1.

    y_0 = np.dot(X, beta)

    noise_levels = np.array([.01, .1, .5, 1., 2., 5., 10., 100., 1000.])
    noise_proto = rng.randn(n_samples)

    noise = noise_levels[np.newaxis, :] * noise_proto[:, np.newaxis]

    y = y_0[:, np.newaxis] + noise
    y = np.hstack([y, noise_proto[:, np.newaxis]])

    ridge_cv.fit(X, y)


