"""Test truncated SVD transformer."""

import numpy as np
import scipy.sparse as sp

import pytest

from sklearn.decomposition import TruncatedSVD, PCA
from sklearn.utils import check_random_state
from sklearn.utils.testing import (assert_array_almost_equal, assert_equal,
                                   assert_raises, assert_greater,
                                   assert_array_less, assert_allclose)


# Make an X that looks somewhat like a small tf-idf matrix.
# XXX newer versions of SciPy have scipy.sparse.rand for this.
shape = 60, 55
n_samples, n_features = shape
rng = check_random_state(42)
X = rng.randint(-100, 20, np.product(shape)).reshape(shape)
X = sp.csr_matrix(np.maximum(X, 0), dtype=np.float64)
X.data[:] = 1 + np.log(X.data)
Xdense = X.A


def test_algorithms():
    svd_a = TruncatedSVD(30, algorithm="arpack")
    svd_r = TruncatedSVD(30, algorithm="randomized", random_state=42)

    Xa = svd_a.fit_transform(X)[:, :6]
    Xr = svd_r.fit_transform(X)[:, :6]
    assert_array_almost_equal(Xa, Xr, decimal=5)

    comp_a = np.abs(svd_a.components_)
    comp_r = np.abs(svd_r.components_)
    # All elements are equal, but some elements are more equal than others.
    assert_array_almost_equal(comp_a[:9], comp_r[:9])
    assert_array_almost_equal(comp_a[9:], comp_r[9:], decimal=2)


def test_attributes():
    for n_components in (10, 25, 41):
        tsvd = TruncatedSVD(n_components).fit(X)
        assert_equal(tsvd.n_components, n_components)
        assert_equal(tsvd.components_.shape, (n_components, n_features))


@pytest.mark.parametrize('algorithm', ("arpack", "randomized"))
def test_too_many_components(algorithm):
    for n_components in (n_features, n_features + 1):
        tsvd = TruncatedSVD(n_components=n_components, algorithm=algorithm)
        assert_raises(ValueError, tsvd.fit, X)


@pytest.mark.parametrize('fmt', ("array", "csr", "csc", "coo", "lil"))
def test_sparse_formats(fmt):
    Xfmt = Xdense if fmt == "dense" else getattr(X, "to" + fmt)()
    tsvd = TruncatedSVD(n_components=11)
    Xtrans = tsvd.fit_transform(Xfmt)
    assert_equal(Xtrans.shape, (n_samples, 11))
    Xtrans = tsvd.transform(Xfmt)
    assert_equal(Xtrans.shape, (n_samples, 11))


@pytest.mark.parametrize('algo', ("arpack", "randomized"))
def test_inverse_transform(algo):
    # We need a lot of components for the reconstruction to be "almost
    # equal" in all positions. XXX Test means or sums instead?
    tsvd = TruncatedSVD(n_components=52, random_state=42, algorithm=algo)
    Xt = tsvd.fit_transform(X)
    Xinv = tsvd.inverse_transform(Xt)
    assert_array_almost_equal(Xinv, Xdense, decimal=1)


def test_integers():
    Xint = X.astype(np.int64)
    tsvd = TruncatedSVD(n_components=6)
    Xtrans = tsvd.fit_transform(Xint)
    assert_equal(Xtrans.shape, (n_samples, tsvd.n_components))


def test_explained_variance():
    # Test sparse data
    svd_a_10_sp = TruncatedSVD(10, algorithm="arpack")
    svd_r_10_sp = TruncatedSVD(10, algorithm="randomized", random_state=42)
    svd_a_20_sp = TruncatedSVD(20, algorithm="arpack")
    svd_r_20_sp = TruncatedSVD(20, algorithm="randomized", random_state=42)
    X_trans_a_10_sp = svd_a_10_sp.fit_transform(X)
    X_trans_r_10_sp = svd_r_10_sp.fit_transform(X)
    X_trans_a_20_sp = svd_a_20_sp.fit_transform(X)
    X_trans_r_20_sp = svd_r_20_sp.fit_transform(X)

    # Test dense data
    svd_a_10_de = TruncatedSVD(10, algorithm="arpack")
    svd_r_10_de = TruncatedSVD(10, algorithm="randomized", random_state=42)
    svd_a_20_de = TruncatedSVD(20, algorithm="arpack")
    svd_r_20_de = TruncatedSVD(20, algorithm="randomized", random_state=42)
    X_trans_a_10_de = svd_a_10_de.fit_transform(X.toarray())
    X_trans_r_10_de = svd_r_10_de.fit_transform(X.toarray())
    X_trans_a_20_de = svd_a_20_de.fit_transform(X.toarray())
    X_trans_r_20_de = svd_r_20_de.fit_transform(X.toarray())

    # helper arrays for tests below
    svds = (svd_a_10_sp, svd_r_10_sp, svd_a_20_sp, svd_r_20_sp, svd_a_10_de,
            svd_r_10_de, svd_a_20_de, svd_r_20_de)
    svds_trans = (
        (svd_a_10_sp, X_trans_a_10_sp),
        (svd_r_10_sp, X_trans_r_10_sp),
        (svd_a_20_sp, X_trans_a_20_sp),
        (svd_r_20_sp, X_trans_r_20_sp),
        (svd_a_10_de, X_trans_a_10_de),
        (svd_r_10_de, X_trans_r_10_de),
        (svd_a_20_de, X_trans_a_20_de),
        (svd_r_20_de, X_trans_r_20_de),
    )
    svds_10_v_20 = (
        (svd_a_10_sp, svd_a_20_sp),
        (svd_r_10_sp, svd_r_20_sp),
        (svd_a_10_de, svd_a_20_de),
        (svd_r_10_de, svd_r_20_de),
    )
    svds_sparse_v_dense = (
        (svd_a_10_sp, svd_a_10_de),
        (svd_a_20_sp, svd_a_20_de),
        (svd_r_10_sp, svd_r_10_de),
        (svd_r_20_sp, svd_r_20_de),
    )

    # Assert the 1st component is equal
    for svd_10, svd_20 in svds_10_v_20:
        assert_array_almost_equal(
            svd_10.explained_variance_ratio_,
            svd_20.explained_variance_ratio_[:10],
            decimal=5,
        )

    # Assert that 20 components has higher explained variance than 10
    for svd_10, svd_20 in svds_10_v_20:
        assert_greater(
            svd_20.explained_variance_ratio_.sum(),
            svd_10.explained_variance_ratio_.sum(),
        )

    # Assert that all the values are greater than 0
    for svd in svds:
        assert_array_less(0.0, svd.explained_variance_ratio_)

    # Assert that total explained variance is less than 1
    for svd in svds:
        assert_array_less(svd.explained_variance_ratio_.sum(), 1.0)

    # Compare sparse vs. dense
    for svd_sparse, svd_dense in svds_sparse_v_dense:
        assert_array_almost_equal(svd_sparse.explained_variance_ratio_,
                                  svd_dense.explained_variance_ratio_)

    # Test that explained_variance is correct
    for svd, transformed in svds_trans:
        total_variance = np.var(X.toarray(), axis=0).sum()
        variances = np.var(transformed, axis=0)
        true_explained_variance_ratio = variances / total_variance

        assert_array_almost_equal(
            svd.explained_variance_ratio_,
            true_explained_variance_ratio,
        )


def test_singular_values():
    # Check that the TruncatedSVD output has the correct singular values

    rng = np.random.RandomState(0)
    n_samples = 100
    n_features = 80

    X = rng.randn(n_samples, n_features)

    apca = TruncatedSVD(n_components=2, algorithm='arpack',
                        random_state=rng).fit(X)
    rpca = TruncatedSVD(n_components=2, algorithm='arpack',
                        random_state=rng).fit(X)
    assert_array_almost_equal(apca.singular_values_, rpca.singular_values_, 12)

    # Compare to the Frobenius norm
    X_apca = apca.transform(X)
    X_rpca = rpca.transform(X)
    assert_array_almost_equal(np.sum(apca.singular_values_**2.0),
                              np.linalg.norm(X_apca, "fro")**2.0, 12)
    assert_array_almost_equal(np.sum(rpca.singular_values_**2.0),
                              np.linalg.norm(X_rpca, "fro")**2.0, 12)

    # Compare to the 2-norms of the score vectors
    assert_array_almost_equal(apca.singular_values_,
                              np.sqrt(np.sum(X_apca**2.0, axis=0)), 12)
    assert_array_almost_equal(rpca.singular_values_,
                              np.sqrt(np.sum(X_rpca**2.0, axis=0)), 12)

    # Set the singular values and see what we get back
    rng = np.random.RandomState(0)
    n_samples = 100
    n_features = 110

    X = rng.randn(n_samples, n_features)

    apca = TruncatedSVD(n_components=3, algorithm='arpack',
                        random_state=rng)
    rpca = TruncatedSVD(n_components=3, algorithm='randomized',
                        random_state=rng)
    X_apca = apca.fit_transform(X)
    X_rpca = rpca.fit_transform(X)

    X_apca /= np.sqrt(np.sum(X_apca**2.0, axis=0))
    X_rpca /= np.sqrt(np.sum(X_rpca**2.0, axis=0))
    X_apca[:, 0] *= 3.142
    X_apca[:, 1] *= 2.718
    X_rpca[:, 0] *= 3.142
    X_rpca[:, 1] *= 2.718

    X_hat_apca = np.dot(X_apca, apca.components_)
    X_hat_rpca = np.dot(X_rpca, rpca.components_)
    apca.fit(X_hat_apca)
    rpca.fit(X_hat_rpca)
    assert_array_almost_equal(apca.singular_values_, [3.142, 2.718, 1.0], 14)
    assert_array_almost_equal(rpca.singular_values_, [3.142, 2.718, 1.0], 14)


def test_truncated_svd_eq_pca():
    # TruncatedSVD should be equal to PCA on centered data

    X_c = X - X.mean(axis=0)

    params = dict(n_components=10, random_state=42)

    svd = TruncatedSVD(algorithm='arpack', **params)
    pca = PCA(svd_solver='arpack', **params)

    Xt_svd = svd.fit_transform(X_c)
    Xt_pca = pca.fit_transform(X_c)

    assert_allclose(Xt_svd, Xt_pca, rtol=1e-9)
    assert_allclose(pca.mean_, 0, atol=1e-9)
    assert_allclose(svd.components_, pca.components_)
