# Author: Tom Dupre la Tour
#
# License: BSD 3 clause

import numpy as np
from numpy.testing import assert_array_equal
from sklearn.utils.testing import assert_allclose
import scipy.sparse as sp

from sklearn.utils.seq_dataset import ArrayDataset64
from sklearn.utils.seq_dataset import CSRDataset64
from sklearn.utils.seq_dataset import ArrayDataset32

from sklearn.datasets import load_iris

iris = load_iris()
X = iris.data.astype(np.float64)
y = iris.target.astype(np.float64)
X_csr = sp.csr_matrix(X)
sample_weight = np.arange(y.size, dtype=np.float64)

X32 = iris.data.astype(np.float32)
y32 = iris.target.astype(np.float32)
X_csr32 = sp.csr_matrix(X32)
sample_weight32 = np.arange(y32.size, dtype=np.float32)


def assert_csr_equal(X, Y):
    X.eliminate_zeros()
    Y.eliminate_zeros()
    assert X.shape[0] == Y.shape[0]
    assert X.shape[1] == Y.shape[1]
    assert_array_equal(X.data, Y.data)
    assert_array_equal(X.indices, Y.indices)
    assert_array_equal(X.indptr, Y.indptr)


def test_seq_dataset():
    dataset1 = ArrayDataset64(X, y, sample_weight, seed=42)
    dataset2 = CSRDataset64(X_csr.data, X_csr.indptr, X_csr.indices,
                          y, sample_weight, seed=42)

    for dataset in (dataset1, dataset2):
        for i in range(5):
            # next sample
            xi_, yi, swi, idx = dataset._next_py()
            xi = sp.csr_matrix((xi_), shape=(1, X.shape[1]))

            assert_csr_equal(xi, X_csr[idx])
            assert yi == y[idx]
            assert swi == sample_weight[idx]

            # random sample
            xi_, yi, swi, idx = dataset._random_py()
            xi = sp.csr_matrix((xi_), shape=(1, X.shape[1]))

            assert_csr_equal(xi, X_csr[idx])
            assert yi == y[idx]
            assert swi == sample_weight[idx]


def test_seq_dataset_shuffle():
    dataset1 = ArrayDataset64(X, y, sample_weight, seed=42)
    dataset2 = CSRDataset64(X_csr.data, X_csr.indptr, X_csr.indices,
                          y, sample_weight, seed=42)

    # not shuffled
    for i in range(5):
        _, _, _, idx1 = dataset1._next_py()
        _, _, _, idx2 = dataset2._next_py()
        assert idx1 == i
        assert idx2 == i

    for i in range(5):
        _, _, _, idx1 = dataset1._random_py()
        _, _, _, idx2 = dataset2._random_py()
        assert idx1 == idx2

    seed = 77
    dataset1._shuffle_py(seed)
    dataset2._shuffle_py(seed)

    for i in range(5):
        _, _, _, idx1 = dataset1._next_py()
        _, _, _, idx2 = dataset2._next_py()
        assert idx1 == idx2

        _, _, _, idx1 = dataset1._random_py()
        _, _, _, idx2 = dataset2._random_py()
        assert idx1 == idx2


def test_fused_types_consistency():
    dataset32 = ArrayDataset32(X32, y32, sample_weight32, seed=42)
    dataset64 = ArrayDataset64(X, y, sample_weight, seed=42)

    for i in range(5):
        # next sample
        xi32, yi32, _, _ = dataset32._next_py()
        xi64, yi64, _, _ = dataset64._next_py()

        xi_data32, _, _ = xi32
        xi_data64, _, _ = xi64

        assert xi_data32.dtype == np.float32
        assert xi_data64.dtype == np.float64
        assert isinstance(yi32, float)
        assert isinstance(yi64, float)

        assert_allclose(xi_data64, xi_data32, rtol=1e-5)
        assert_allclose(yi64, yi32, rtol=1e-5)
