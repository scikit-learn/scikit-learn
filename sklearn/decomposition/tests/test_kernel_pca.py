from datetime import datetime

import numpy as np
import scipy.sparse as sp
import pytest
from sklearn.base import BaseEstimator, TransformerMixin

from sklearn.exceptions import KernelWarning
from sklearn.utils.validation import check_kernel_eigenvalues
from sklearn.utils.testing import (assert_array_almost_equal, assert_less,
                                   assert_equal, assert_not_equal,
                                   assert_raises, assert_warns,
                                   assert_no_warnings)

from sklearn.decomposition import PCA, KernelPCA
from sklearn.datasets import make_circles
from sklearn.linear_model import Perceptron
from sklearn.pipeline import Pipeline
from sklearn.model_selection import GridSearchCV
from sklearn.metrics.pairwise import rbf_kernel


def test_kernel_pca():
    """ Nominal test for all solvers and all known kernels + a custom one.
    It tests
     - that fit_transform is equivalent to fit+transform
     - that the shapes of transforms and inverse transforms are correct """
    rng = np.random.RandomState(0)
    X_fit = rng.random_sample((5, 4))
    X_pred = rng.random_sample((2, 4))

    def histogram(x, y, **kwargs):
        # Histogram kernel implemented as a callable.
        assert_equal(kwargs, {})    # no kernel_params that we didn't ask for
        return np.minimum(x, y).sum()

    for eigen_solver in ("auto", "dense", "arpack", "randomized"):
        for kernel in ("linear", "rbf", "poly", histogram):
            # histogram kernel produces singular matrix inside linalg.solve
            # XXX use a least-squares approximation?
            inv = not callable(kernel)

            # transform fit data
            kpca = KernelPCA(4, kernel=kernel, eigen_solver=eigen_solver,
                             fit_inverse_transform=inv)
            X_fit_transformed = kpca.fit_transform(X_fit)
            X_fit_transformed2 = kpca.fit(X_fit).transform(X_fit)
            assert_array_almost_equal(np.abs(X_fit_transformed),
                                      np.abs(X_fit_transformed2))

            # non-regression test: previously, gamma would be 0 by default,
            # forcing all eigenvalues to 0 under the poly kernel
            assert_not_equal(X_fit_transformed.size, 0)

            # transform new data
            X_pred_transformed = kpca.transform(X_pred)
            assert_equal(X_pred_transformed.shape[1],
                         X_fit_transformed.shape[1])

            # inverse transform
            if inv:
                X_pred2 = kpca.inverse_transform(X_pred_transformed)
                assert_equal(X_pred2.shape, X_pred.shape)


def test_kernel_pca_invalid_parameters():
    assert_raises(ValueError, KernelPCA, 10, fit_inverse_transform=True,
                  kernel='precomputed')


def test_kernel_pca_consistent_transform():
    """ Tests that after fitting a kPCA model, it is independent of the
    original data object (uses an inner copy) """
    # X_fit_ needs to retain the old, unmodified copy of X
    state = np.random.RandomState(0)
    X = state.rand(10, 10)
    kpca = KernelPCA(random_state=state).fit(X)
    transformed1 = kpca.transform(X)

    X_copy = X.copy()
    X[:, 0] = 666
    transformed2 = kpca.transform(X_copy)
    assert_array_almost_equal(transformed1, transformed2)


def test_kernel_pca_sparse():
    """ Tests that kPCA works on a sparse data input. Same test than
    test_kernel_pca except inverse_transform (why?) """
    rng = np.random.RandomState(0)
    X_fit = sp.csr_matrix(rng.random_sample((5, 4)))
    X_pred = sp.csr_matrix(rng.random_sample((2, 4)))

    for eigen_solver in ("auto", "arpack", "randomized"):
        for kernel in ("linear", "rbf", "poly"):
            # transform fit data
            kpca = KernelPCA(4, kernel=kernel, eigen_solver=eigen_solver,
                             fit_inverse_transform=False)
            X_fit_transformed = kpca.fit_transform(X_fit)
            X_fit_transformed2 = kpca.fit(X_fit).transform(X_fit)
            assert_array_almost_equal(np.abs(X_fit_transformed),
                                      np.abs(X_fit_transformed2))

            # transform new data
            X_pred_transformed = kpca.transform(X_pred)
            assert_equal(X_pred_transformed.shape[1],
                         X_fit_transformed.shape[1])

            # inverse transform
            # X_pred2 = kpca.inverse_transform(X_pred_transformed)
            # assert_equal(X_pred2.shape, X_pred.shape)


def test_kernel_pca_linear_kernel():
    """ Tests that kPCA with a linear kernel is equivalent to PCA """
    rng = np.random.RandomState(0)
    X_fit = rng.random_sample((5, 4))
    X_pred = rng.random_sample((2, 4))

    # for a linear kernel, kernel PCA should find the same projection as PCA
    # modulo the sign (direction)
    # fit only the first four components: fifth is near zero eigenvalue, so
    # can be trimmed due to roundoff error
    assert_array_almost_equal(
        np.abs(KernelPCA(4).fit(X_fit).transform(X_pred)),
        np.abs(PCA(4).fit(X_fit).transform(X_pred)))


def test_kernel_pca_linear_kernel2():
    """ Tests that kPCA with a linear kernel is equivalent to PCA, for all
    solvers"""
    rng = np.random.RandomState(0)
    X_fit = rng.random_sample((6, 10))
    X_pred = rng.random_sample((2, 10))

    # for a linear kernel, kernel PCA should find the same projection as PCA
    # modulo the sign (direction)
    for solver in ("auto", "dense", "arpack", "randomized"):
        assert_array_almost_equal(
            np.abs(KernelPCA(4, eigen_solver=solver).fit(X_fit)
                   .transform(X_pred)),
            np.abs(PCA(4, svd_solver=solver if solver != "dense" else "full")
                   .fit(X_fit).transform(X_pred)))


def test_kernel_pca_n_components():
    """ Tests that the number of components selected is correctly taken into
    account for projections, for all solvers """
    rng = np.random.RandomState(0)
    X_fit = rng.random_sample((5, 4))
    X_pred = rng.random_sample((2, 4))

    for eigen_solver in ("dense", "arpack", "randomized"):
        for c in [1, 2, 4]:
            kpca = KernelPCA(n_components=c, eigen_solver=eigen_solver)
            shape = kpca.fit(X_fit).transform(X_pred).shape

            assert_equal(shape, (2, c))


def test_remove_zero_eig():
    """ Tests that the null-space (Zero) eigenvalues are removed when
    remove_zero_eig=True, whereas they are not by default """
    X = np.array([[1 - 1e-30, 1], [1, 1], [1, 1 - 1e-20]])

    # n_components=None (default) => remove_zero_eig is True
    kpca = KernelPCA()
    Xt = kpca.fit_transform(X)
    assert_equal(Xt.shape, (3, 0))

    kpca = KernelPCA(n_components=2)
    Xt = kpca.fit_transform(X)
    assert_equal(Xt.shape, (3, 2))

    kpca = KernelPCA(n_components=2, remove_zero_eig=True)
    Xt = kpca.fit_transform(X)
    assert_equal(Xt.shape, (3, 0))


def test_leave_zero_eig():
    """This test checks that fit().transform() returns the same result than
    fit_transform() in case of non-removed zero eigenvalue"""
    X_fit = np.array([[1, 1], [0, 0]])

    # Assert that even with all np warnings on, there is no div by zero warning
    with pytest.warns(None) as record:
        with np.errstate(all='warn'):
            k = KernelPCA(n_components=2, remove_zero_eig=False,
                          eigen_solver="dense")
            # Fit, then transform
            A = k.fit(X_fit).transform(X_fit)
            # Do both at once
            B = k.fit_transform(X_fit)
            # Compare
            assert_array_almost_equal(np.abs(A), np.abs(B))

    for w in record:
        # There might be warnings about the kernel being badly conditioned,
        # but there should not be warnings about division by zero.
        # (Numpy division by zero warning can have many message variants, but
        # at least we know that it is a RuntimeWarning so lets check only this)
        assert w.category is not RuntimeWarning


def test_kernel_pca_precomputed():
    """ Tests that kPCA works when the kernel has been precomputed, for all
    solvers """
    rng = np.random.RandomState(0)
    X_fit = rng.random_sample((5, 4))
    X_pred = rng.random_sample((2, 4))

    for eigen_solver in ("dense", "arpack", "randomized"):
        X_kpca = KernelPCA(4, eigen_solver=eigen_solver).\
            fit(X_fit).transform(X_pred)
        X_kpca2 = KernelPCA(
            4, eigen_solver=eigen_solver, kernel='precomputed').fit(
                np.dot(X_fit, X_fit.T)).transform(np.dot(X_pred, X_fit.T))

        X_kpca_train = KernelPCA(
            4, eigen_solver=eigen_solver,
            kernel='precomputed').fit_transform(np.dot(X_fit, X_fit.T))
        X_kpca_train2 = KernelPCA(
            4, eigen_solver=eigen_solver, kernel='precomputed').fit(
                np.dot(X_fit, X_fit.T)).transform(np.dot(X_fit, X_fit.T))

        assert_array_almost_equal(np.abs(X_kpca),
                                  np.abs(X_kpca2))

        assert_array_almost_equal(np.abs(X_kpca_train),
                                  np.abs(X_kpca_train2))


def test_kernel_pca_invalid_kernel():
    """ Tests that using an invalid kernel name raises a ValueError at fit
    time"""
    rng = np.random.RandomState(0)
    X_fit = rng.random_sample((2, 4))
    kpca = KernelPCA(kernel="tototiti")
    assert_raises(ValueError, kpca.fit, X_fit)


@pytest.mark.filterwarnings('ignore: The default of the `iid`')  # 0.22
def test_gridsearch_pipeline():
    # Test if we can do a grid-search to find parameters to separate
    # circles with a perceptron model.
    X, y = make_circles(n_samples=400, factor=.3, noise=.05,
                        random_state=0)
    kpca = KernelPCA(kernel="rbf", n_components=2)
    pipeline = Pipeline([("kernel_pca", kpca),
                         ("Perceptron", Perceptron(max_iter=5))])
    param_grid = dict(kernel_pca__gamma=2. ** np.arange(-2, 2))
    grid_search = GridSearchCV(pipeline, cv=3, param_grid=param_grid)
    grid_search.fit(X, y)
    assert_equal(grid_search.best_score_, 1)


@pytest.mark.filterwarnings('ignore: The default of the `iid`')  # 0.22
def test_gridsearch_pipeline_precomputed():
    # Test if we can do a grid-search to find parameters to separate
    # circles with a perceptron model using a precomputed kernel.
    X, y = make_circles(n_samples=400, factor=.3, noise=.05,
                        random_state=0)
    kpca = KernelPCA(kernel="precomputed", n_components=2)
    pipeline = Pipeline([("kernel_pca", kpca),
                         ("Perceptron", Perceptron(max_iter=5))])
    param_grid = dict(Perceptron__max_iter=np.arange(1, 5))
    grid_search = GridSearchCV(pipeline, cv=3, param_grid=param_grid)
    X_kernel = rbf_kernel(X, gamma=2.)
    grid_search.fit(X_kernel, y)
    assert_equal(grid_search.best_score_, 1)


def test_nested_circles():
    # Test the linear separability of the first 2D KPCA transform
    X, y = make_circles(n_samples=400, factor=.3, noise=.05,
                        random_state=0)

    # 2D nested circles are not linearly separable
    train_score = Perceptron(max_iter=5).fit(X, y).score(X, y)
    assert_less(train_score, 0.8)

    # Project the circles data into the first 2 components of a RBF Kernel
    # PCA model.
    # Note that the gamma value is data dependent. If this test breaks
    # and the gamma value has to be updated, the Kernel PCA example will
    # have to be updated too.
    kpca = KernelPCA(kernel="rbf", n_components=2,
                     fit_inverse_transform=True, gamma=2.)
    X_kpca = kpca.fit_transform(X)

    # The data is perfectly linearly separable in that space
    train_score = Perceptron(max_iter=5).fit(X_kpca, y).score(X_kpca, y)
    assert_equal(train_score, 1.0)


def test_errors_and_warnings():
    """Tests that bad kernels raise error and warnings"""

    solvers = ['dense', 'arpack', 'randomized']
    solvers_except_arpack = ['dense', 'randomized']

    # First create an identity transformer class
    # ------------------------------------------
    class IdentityKernelTransformer(BaseEstimator, TransformerMixin):
        """We will use this transformer so that the passed kernel matrix is
        not centered when kPCA is fit"""

        def __init__(self):
            pass

        def fit(self, K, y=None):
            return self

        def transform(self, K, y=None, copy=True):
            return K

    # Significant imaginary parts: error
    # ----------------------------------
    # As of today it seems that the kernel matrix is always cast as a float
    # whatever the method (precomputed or callable).
    # The following test therefore fails ("did not raise").
    K = [[5, 0],
         [0, 6 * 1e-5j]]
    # for solver in solvers:
    #     kpca = KernelPCA(kernel=kernel_getter, eigen_solver=solver,
    #                      fit_inverse_transform=False)
    #     kpca._centerer = IdentityKernelTransformer()
    #     K = kpca._get_kernel(K)
    #     with pytest.raises(ValueError):
    #         # note: we can not test 'fit' because _centerer would be replaced
    #         kpca._fit_transform(K)
    #
    # For safety concerning future evolutions the corresponding code is left in
    # KernelPCA, and we test it directly by calling the inner method here:
    with pytest.raises(ValueError):
        check_kernel_eigenvalues((K[0][0], K[1][1]))

    # All negative eigenvalues: error
    # -------------------------------
    K = [[-5, 0],
         [0, -6e-5]]
    # check that the inner method works
    with pytest.raises(ValueError):
        check_kernel_eigenvalues((K[0][0], K[1][1]))

    for solver in solvers:
        kpca = KernelPCA(kernel="precomputed", eigen_solver=solver,
                         fit_inverse_transform=False)
        kpca._centerer = IdentityKernelTransformer()
        K = kpca._get_kernel(K)
        with pytest.raises(ValueError):
            # note: we can not test 'fit' because _centerer would be replaced
            kpca._fit_transform(K)

    # Significant negative eigenvalue: warning
    # ----------------------------------------
    K = [[5, 0],
         [0, -6e-5]]
    # check that the inner method works
    assert_warns(KernelWarning,
                 lambda: check_kernel_eigenvalues((K[0][0], K[1][1])))

    for solver in solvers_except_arpack:
        # Note: arpack detects this case and raises an error already
        kpca = KernelPCA(kernel="precomputed", eigen_solver=solver,
                         fit_inverse_transform=False)
        kpca._centerer = IdentityKernelTransformer()
        K = kpca._get_kernel(K)
        # note: we can not test 'fit' because _centerer would be replaced
        assert_warns(KernelWarning, lambda: kpca._fit_transform(K))

    # Bad conditioning
    # ----------------
    K = [[5, 0],
         [0, 4e-12]]
    # check that the inner method works
    assert_warns(KernelWarning,
                 lambda: check_kernel_eigenvalues((K[0][0], K[1][1]),
                                                  warn_on_zeros=True))

    # This is actually a normal behaviour when the number of samples is big,
    # so no special warning should be raised in this case
    for solver in solvers_except_arpack:
        # Note: arpack detects this case and raises an error already
        kpca = KernelPCA(kernel="precomputed", eigen_solver=solver,
                         fit_inverse_transform=False)
        kpca._centerer = IdentityKernelTransformer()
        K = kpca._get_kernel(K)
        # note: we can not test 'fit' because _centerer would be replaced
        assert_no_warnings(lambda: kpca._fit_transform(K))


def test_kernel_pca_time_and_equivalence():
    """Checks that 'dense', 'arpack' and 'randomized' solvers give similar
    results and benchmarks their respective execution times. This test can
    be transformed into a benchmark by setting benchmark_mode to True.
    """

    # Generate random data
    n_training_samples = 2000
    n_features = 10
    rng = np.random.RandomState(0)
    X_fit = rng.random_sample((n_training_samples, n_features))
    X_pred = rng.random_sample((100, n_features))

    # Experiments design
    n_compo_range = [4, 10, 20]
    arpack_all = False
    n_iter = 3

    # Runs
    ref_time = np.empty((len(n_compo_range), n_iter)) * np.nan
    a_time = np.empty((len(n_compo_range), n_iter)) * np.nan
    r_time = np.empty((len(n_compo_range), n_iter)) * np.nan
    for j, n_components in enumerate(n_compo_range):

        # reference (full)
        for i in range(n_iter):
            start_time = datetime.now()
            ref_pred = KernelPCA(n_components, eigen_solver="dense")\
                .fit(X_fit).transform(X_pred)
            ref_time[j, i] = (datetime.now() - start_time).total_seconds()

        # arpack
        if arpack_all or n_components < 100:
            for i in range(n_iter):
                start_time = datetime.now()
                a_pred = KernelPCA(n_components, eigen_solver="arpack")\
                    .fit(X_fit).transform(X_pred)
                # check that the result is still correct despite the approx
                assert_array_almost_equal(np.abs(a_pred), np.abs(ref_pred))
                a_time[j, i] = (datetime.now() - start_time).total_seconds()

        # randomized
        for i in range(n_iter):
            start_time = datetime.now()
            r_pred = KernelPCA(n_components, eigen_solver="randomized")\
                .fit(X_fit).transform(X_pred)
            # check that the result is still correct despite the approximation
            assert_array_almost_equal(np.abs(r_pred), np.abs(ref_pred))
            r_time[j, i] = (datetime.now() - start_time).total_seconds()

    # Compute statistics for the 3 methods
    avg_ref_time = ref_time.mean(axis=1)
    avg_a_time = a_time.mean(axis=1)
    avg_r_time = r_time.mean(axis=1)

    # A few asserts - very "light" because Continuous integration target
    # platforms vary quite a lot.
    #
    # Check that randomized reduces the time on average
    assert np.mean(avg_r_time / avg_ref_time) < 0.75

    # Check that arpack reduces the time at least on one run
    assert min(avg_a_time / avg_ref_time) < 0.75
