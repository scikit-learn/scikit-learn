"""Test latent semantic analysis."""

import numpy as np
import scipy.sparse as sp

from sklearn.decomposition import LSA
from sklearn.utils import check_random_state
from sklearn.utils.testing import assert_equal, assert_raises


# Make an X that looks somewhat like a small tf-idf matrix.
# XXX newer versions of SciPy have scipy.sparse.rand for this.
shape = 60, 55
n_samples, n_features = shape
rng = check_random_state(42)
X = rng.randint(-100, 20, np.product(shape)).reshape(shape)
X = sp.csr_matrix(np.maximum(X, 0), dtype=np.float64)
X.data[:] = 1 + np.log(X.data)


def test_attributes():
    for n_components in (10, 25, 41):
        lsa = LSA(n_components).fit(X)
        assert_equal(lsa.n_components, n_components)
        assert_equal(lsa.components_.shape, (n_components, n_features))


def test_too_many_components():
    for n_components in (n_features, n_features+1):
        lsa = LSA(n_components=n_components)
        assert_raises(ValueError, lsa.fit, X)


def test_sparse_formats():
    for fmt in ("array", "csr", "csc", "coo", "lil"):
        Xfmt = getattr(X, "to" + fmt)()
        lsa = LSA(n_components=11)
        Xtrans = lsa.fit_transform(Xfmt)
        assert_equal(Xtrans.shape, (n_samples, 11))
        Xtrans = lsa.transform(Xfmt)
        assert_equal(Xtrans.shape, (n_samples, 11))


def test_integers():
    Xint = X.astype(np.int64)
    lsa = LSA(n_components=6)
    Xtrans = lsa.fit_transform(Xint)
    assert_equal(Xtrans.shape, (n_samples, lsa.n_components))
