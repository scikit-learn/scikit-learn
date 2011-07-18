import numpy as np

from numpy.testing import assert_array_almost_equal
from nose.tools import assert_equal
from nose.tools import assert_raises
from nose import SkipTest

from .. import PCA, KernelPCA


def test_kernel_pca():
    X_fit = np.random.random((5, 4))
    X_pred = np.random.random((2, 4))

    for kernel in ("linear", "rbf", "poly"):
        # transform fit data
        kpca = KernelPCA(kernel=kernel, fit_inverse_transform=True)
        X_fit_transformed = kpca.fit_transform(X_fit)
        X_fit_transformed2 = kpca.fit(X_fit).transform(X_fit)
        assert_array_almost_equal(X_fit_transformed, X_fit_transformed2)

        # transform new data
        X_pred_transformed = kpca.transform(X_pred)
        assert_equal(X_pred_transformed.shape[1], X_fit_transformed.shape[1])

        # inverse transform
        X_pred2 = kpca.inverse_transform(X_pred_transformed)
        assert_equal(X_pred2.shape, X_pred.shape)


def test_kernel_pca_linear_kernel():
    raise SkipTest
    X_fit = np.random.random((5, 4))
    X_pred = np.random.random((2, 4))

    # for a linear kernel, kernel PCA should find the same projection as PCA
    # modulo the sign (direction)
    assert_array_almost_equal(np.abs(KernelPCA().fit(X_fit).transform(X_pred)),
                              np.abs(PCA().fit(X_fit).transform(X_pred)))


def test_kernel_pca_n_components():
    X_fit = np.random.random((5, 4))
    X_pred = np.random.random((2, 4))

    for c in [1, 2, 4]:
        kpca = KernelPCA(n_components = c)
        shape = kpca.fit(X_fit).transform(X_pred).shape

        assert_equal(shape, (2, c))


def test_kernel_pca_precomputed():
    X_fit = np.random.random((5, 4))
    X_pred = np.random.random((2, 4))

    X_kpca = KernelPCA().fit(X_fit).transform(X_pred)
    X_kpca2 = KernelPCA(kernel="precomputed").fit(np.dot(X_fit, X_fit.T)). \
              transform(np.dot(X_pred, X_fit.T))

    assert_array_almost_equal(X_kpca, X_kpca2)


def test_kernel_pca_invalid_kernel():
    X_fit = np.random.random((2, 4))
    kpca = KernelPCA(kernel="tototiti")
    assert_raises(ValueError, kpca.fit, X_fit)