"""Test truncated SVD transformer."""

import numpy as np
import scipy.sparse as sp

from sklearn.decomposition import TruncatedSVD
from sklearn.utils import check_random_state
from sklearn.utils.testing import (assert_array_almost_equal, assert_equal,
                                   assert_raises)


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
    assert_array_almost_equal(Xa, Xr)

    comp_a = np.abs(svd_a.components_)
    comp_r = np.abs(svd_r.components_)
    # All elements are equal, but some elements are more equal than others.
    assert_array_almost_equal(comp_a[:9], comp_r[:9])
    assert_array_almost_equal(comp_a[9:], comp_r[9:], decimal=3)


def test_attributes():
    for n_components in (10, 25, 41):
        tsvd = TruncatedSVD(n_components).fit(X)
        assert_equal(tsvd.n_components, n_components)
        assert_equal(tsvd.components_.shape, (n_components, n_features))


def test_too_many_components():
    for algorithm in ["arpack", "randomized"]:
        for n_components in (n_features, n_features+1):
            tsvd = TruncatedSVD(n_components=n_components, algorithm=algorithm)
            assert_raises(ValueError, tsvd.fit, X)


def test_sparse_formats():
    for fmt in ("array", "csr", "csc", "coo", "lil"):
        Xfmt = Xdense if fmt == "dense" else getattr(X, "to" + fmt)()
        tsvd = TruncatedSVD(n_components=11)
        Xtrans = tsvd.fit_transform(Xfmt)
        assert_equal(Xtrans.shape, (n_samples, 11))
        Xtrans = tsvd.transform(Xfmt)
        assert_equal(Xtrans.shape, (n_samples, 11))


def test_inverse_transform():
    for algo in ("arpack", "randomized"):
        # We need a lot of components for the reconstruction to be "almost
        # equal" in all positions. XXX Test means or sums instead?
        tsvd = TruncatedSVD(n_components=52, random_state=42)
        Xt = tsvd.fit_transform(X)
        Xinv = tsvd.inverse_transform(Xt)
        assert_array_almost_equal(Xinv, Xdense, decimal=1)


def test_integers():
    Xint = X.astype(np.int64)
    tsvd = TruncatedSVD(n_components=6)
    Xtrans = tsvd.fit_transform(Xint)
    assert_equal(Xtrans.shape, (n_samples, tsvd.n_components))
