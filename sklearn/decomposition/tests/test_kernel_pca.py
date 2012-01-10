import numpy as np
import scipy.sparse as sp

from numpy.testing import assert_array_almost_equal
from nose.tools import assert_equal
from nose.tools import assert_raises

from .. import PCA, KernelPCA


def test_kernel_pca():
    rng = np.random.RandomState(0)
    X_fit = rng.random_sample((5, 4))
    X_pred = rng.random_sample((2, 4))

    for eigen_solver in ("auto", "dense", "arpack"):
        for kernel in ("linear", "rbf", "poly"):
            # transform fit data
            kpca = KernelPCA(4, kernel=kernel, eigen_solver=eigen_solver,
                             fit_inverse_transform=True)
            X_fit_transformed = kpca.fit_transform(X_fit)
            X_fit_transformed2 = kpca.fit(X_fit).transform(X_fit)
            assert_array_almost_equal(np.abs(X_fit_transformed),
                                      np.abs(X_fit_transformed2))

            # transform new data
            X_pred_transformed = kpca.transform(X_pred)
            assert_equal(X_pred_transformed.shape[1],
                         X_fit_transformed.shape[1])

            # inverse transform
            X_pred2 = kpca.inverse_transform(X_pred_transformed)
            assert_equal(X_pred2.shape, X_pred.shape)


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
            #X_pred2 = kpca.inverse_transform(X_pred_transformed)
            #assert_equal(X_pred2.shape, X_pred.shape)


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


def test_kernel_pca_precomputed():
    rng = np.random.RandomState(0)
    X_fit = rng.random_sample((5, 4))
    X_pred = rng.random_sample((2, 4))

    for eigen_solver in ("dense", "arpack"):
        X_kpca = KernelPCA(4, eigen_solver=eigen_solver).\
            fit(X_fit).transform(X_pred)
        X_kpca2 = KernelPCA(4, kernel="precomputed",
                            eigen_solver=eigen_solver).\
                            fit(np.dot(X_fit, X_fit.T)).\
                            transform(np.dot(X_pred, X_fit.T))

        assert_array_almost_equal(np.abs(X_kpca),
                                  np.abs(X_kpca2))


def test_kernel_pca_invalid_kernel():
    rng = np.random.RandomState(0)
    X_fit = rng.random_sample((2, 4))
    kpca = KernelPCA(kernel="tototiti")
    assert_raises(ValueError, kpca.fit, X_fit)


if __name__ == '__main__':
    import nose
    nose.run(argv=['', __file__])
