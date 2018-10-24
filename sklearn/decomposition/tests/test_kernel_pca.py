import numpy as np
import scipy.sparse as sp
import pytest

from sklearn.utils.testing import (assert_array_almost_equal, assert_less,
                                   assert_equal, assert_not_equal,
                                   assert_raises, ignore_warnings)

from sklearn.decomposition import PCA, KernelPCA
from sklearn.datasets import make_circles
from sklearn.linear_model import Perceptron
from sklearn.pipeline import Pipeline
from sklearn.model_selection import GridSearchCV
from sklearn.metrics.pairwise import rbf_kernel


def test_kernel_pca():
    rng = np.random.RandomState(0)
    X_fit = rng.random_sample((5, 4))
    X_pred = rng.random_sample((2, 4))

    def histogram(x, y, **kwargs):
        # Histogram kernel implemented as a callable.
        assert_equal(kwargs, {})    # no kernel_params that we didn't ask for
        return np.minimum(x, y).sum()

    for eigen_solver in ("auto", "dense", "arpack"):
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
    rng = np.random.RandomState(0)
    X_fit = sp.csr_matrix(rng.random_sample((5, 4)))
    X_pred = sp.csr_matrix(rng.random_sample((2, 4)))

    for eigen_solver in ("auto", "arpack"):
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


def test_kernel_pca_n_components():
    rng = np.random.RandomState(0)
    X_fit = rng.random_sample((5, 4))
    X_pred = rng.random_sample((2, 4))

    for eigen_solver in ("dense", "arpack"):
        for c in [1, 2, 4]:
            kpca = KernelPCA(n_components=c, eigen_solver=eigen_solver)
            shape = kpca.fit(X_fit).transform(X_pred).shape

            assert_equal(shape, (2, c))


def test_remove_zero_eig():
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


def test_kernel_pca_precomputed():
    rng = np.random.RandomState(0)
    X_fit = rng.random_sample((5, 4))
    X_pred = rng.random_sample((2, 4))

    for eigen_solver in ("dense", "arpack"):
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
    rng = np.random.RandomState(0)
    X_fit = rng.random_sample((2, 4))
    kpca = KernelPCA(kernel="tototiti")
    assert_raises(ValueError, kpca.fit, X_fit)


@pytest.mark.filterwarnings('ignore: The default of the `iid`')  # 0.22
# 0.23. warning about tol not having its correct default value.
@pytest.mark.filterwarnings('ignore:max_iter and tol parameters have been')
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
# 0.23. warning about tol not having its correct default value.
@pytest.mark.filterwarnings('ignore:max_iter and tol parameters have been')
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


# 0.23. warning about tol not having its correct default value.
@pytest.mark.filterwarnings('ignore:max_iter and tol parameters have been')
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
