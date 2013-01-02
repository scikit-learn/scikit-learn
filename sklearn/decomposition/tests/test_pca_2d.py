import numpy as np


from numpy.testing import assert_array_almost_equal
from numpy.testing import assert_array_equal

from sklearn import datasets
from sklearn.decomposition import PCA2D


digits = datasets.load_digits()


def test_correct_shapes():

    rng = np.random.RandomState(0)

    X = rng.randn(50, 30, 40)

    X_copy = np.copy(X)

    pca2d = PCA2D(n_row_components=10, n_column_components=12)

    # Testing shapes after transformation

    pca2d.fit(X)

    assert_array_equal(pca2d.row_components_.shape, (30, 10))

    assert_array_equal(pca2d.column_components_.shape, (40, 12))

    Xr = pca2d.transform(X)

    assert_array_equal(Xr.shape, (50, 10, 12))

    #Xr_2d = pca2d.transform(X[0, :, :])

    #assert_array_equal(Xr_2d.shape, (10, 12))

    # Testing shapes after inverse transformation

    # Testing first that the data X has not been changed by the transformation

    assert_array_almost_equal(X, X_copy)

    X_inv = pca2d.inverse_transform(Xr)

    assert_array_equal(X_inv.shape, X.shape)

    #X_inv_2d = pca2d.inverse_transform(Xr_2d)

    #assert_array_equal(X_inv_2d.shape, X[0, :, :].shape)


def test_logic():

    # Testing that after (only row and column)transformation there is no
    # correlation between the new features
    X = digits.images

    # after row transformation

    pca2d = PCA2D(n_row_components=4, n_column_components=0, epsilon=0)
    pca2d.fit(X)
    Xr = pca2d.transform(X)
    # Computing the covariance matrix of Xr on the row dimension
    X_cov = np.tensordot(Xr, Xr, axes=([0, 2], [0, 2])) / X.shape[0]
    X_cov = X_cov - np.diag(np.diag(X_cov))
    assert_array_almost_equal(X_cov, np.zeros((4, 4)))

    # after column transformation
    pca2d = PCA2D(n_row_components=0, n_column_components=4, epsilon=0)
    pca2d.fit(X)
    Xr = pca2d.transform(X)
    # Computing the covariance matrix of Xr on the column dimension
    X_cov = np.tensordot(Xr, Xr, axes=([0, 1], [0, 1])) / X.shape[0]
    X_cov = X_cov - np.diag(np.diag(X_cov))
    assert_array_almost_equal(X_cov, np.zeros((4, 4)))


def test_logic_whitening():

    # Testing that after transformation and whitening (row and column)
    # we have unit variance for the new features

    X = digits.images
    # after row transformation
    pca2d = PCA2D(n_row_components=4, n_column_components=0,
                  epsilon=0, row_whiten=True)
    pca2d.fit(X)
    Xr = pca2d.transform(X)
    # Computing the covariance matrix of Xr on the row dimension
    X_cov = np.tensordot(Xr, Xr, axes=([0, 2], [0, 2])) / X.shape[0]
    assert_array_almost_equal(X_cov, np.eye(4, 4))

    # after column transformation
    pca2d = PCA2D(n_row_components=0, n_column_components=4,
                  epsilon=0, column_whiten=True)
    pca2d.fit(X)
    Xr = pca2d.transform(X)
    # Computing the covariance matrix of Xr on the column dimension
    X_cov = np.tensordot(Xr, Xr, axes=([0, 1], [0, 1])) / X.shape[0]
    assert_array_almost_equal(X_cov, np.eye(4, 4))


def test_pca2d_efficiency():

    X = np.arange(2000 * 20 * 20).reshape(2000, 20, 20)
    # Testing performance when 90% of the variance is kept (row and column
    pca2d = PCA2D(n_row_components=0.9, n_column_components=0.9)
    pca2d.fit(X)
    Xr = pca2d.transform(X)

    # Due to the high correlation between the features, only one is
    # needed to get 90% of the variance
    assert_array_equal(Xr.shape, (2000, 1, 1))

    X_ori = pca2d.inverse_transform(Xr)

    # X_ori should be approximately equal to X
    assert_array_almost_equal(X_ori, X)
