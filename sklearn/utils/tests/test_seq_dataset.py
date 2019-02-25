# Author: Tom Dupre la Tour
#         Joan Massich <mailsik@gmail.com>
#
# License: BSD 3 clause

import pytest
import numpy as np
from numpy.testing import assert_array_equal
from sklearn.utils.testing import assert_allclose
import scipy.sparse as sp

from sklearn.utils.seq_dataset import ArrayDataset64
from sklearn.utils.seq_dataset import ArrayDataset32
from sklearn.utils.seq_dataset import CSRDataset64
from sklearn.utils.seq_dataset import CSRDataset32

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


def assert_csr_equal_values(current, expected):
    current.eliminate_zeros()
    expected.eliminate_zeros()
    expected = expected.astype(current.dtype)
    assert current.shape[0] == expected.shape[0]
    assert current.shape[1] == expected.shape[1]
    assert_array_equal(current.data, expected.data)
    assert_array_equal(current.indices, expected.indices)
    assert_array_equal(current.indptr, expected.indptr)


@pytest.mark.parametrize(
    'dataset',
    [
        ArrayDataset64(X, y, sample_weight, seed=42),
        ArrayDataset32(X32, y32, sample_weight32, seed=42),
        CSRDataset64(X_csr.data, X_csr.indptr, X_csr.indices, y, sample_weight,
                     seed=42),
        CSRDataset32(X_csr32.data, X_csr32.indptr, X_csr32.indices, y32,
                     sample_weight32, seed=42),
    ],
    ids=['ArrayDataset64', 'ArrayDataset32', 'CSRDataset64', 'CSRDataset32']
)
def test_seq_dataset_basic_iteration(dataset):
    NUMBER_OF_RUNS = 5
    for _ in range(NUMBER_OF_RUNS):
        # next sample
        xi_, yi, swi, idx = dataset._next_py()
        xi = sp.csr_matrix((xi_), shape=(1, X.shape[1]))

        assert_csr_equal_values(xi, X_csr[idx])
        assert yi == y[idx]
        assert swi == sample_weight[idx]

        # random sample
        xi_, yi, swi, idx = dataset._random_py()
        xi = sp.csr_matrix((xi_), shape=(1, X.shape[1]))

        assert_csr_equal_values(xi, X_csr[idx])
        assert yi == y[idx]
        assert swi == sample_weight[idx]


@pytest.mark.parametrize(
    'dense_dataset,sparse_dataset',
    [
        (ArrayDataset64(X, y, sample_weight, seed=42),
        CSRDataset64(X_csr.data, X_csr.indptr, X_csr.indices, y, sample_weight,
                     seed=42)),
        (ArrayDataset32(X32, y32, sample_weight32, seed=42),
        CSRDataset32(X_csr32.data, X_csr32.indptr, X_csr32.indices, y32,
                     sample_weight32, seed=42))
    ],
    ids=['float 64', 'float 32']
)
def test_seq_dataset_shuffle(dense_dataset, sparse_dataset):
    # not shuffled
    for i in range(5):
        _, _, _, idx1 = dense_dataset._next_py()
        _, _, _, idx2 = sparse_dataset._next_py()
        assert idx1 == i
        assert idx2 == i

    for i in range(5):
        _, _, _, idx1 = dense_dataset._random_py()
        _, _, _, idx2 = sparse_dataset._random_py()
        assert idx1 == idx2

    seed = 77
    dense_dataset._shuffle_py(seed)
    sparse_dataset._shuffle_py(seed)

    for i in range(5):
        _, _, _, idx1 = dense_dataset._next_py()
        _, _, _, idx2 = sparse_dataset._next_py()
        assert idx1 == idx2

        _, _, _, idx1 = dense_dataset._random_py()
        _, _, _, idx2 = sparse_dataset._random_py()
        assert idx1 == idx2


@pytest.mark.parametrize(
    'dataset32,dataset64',
    [
        (
            ArrayDataset32(X32, y32, sample_weight32, seed=42),
            ArrayDataset64(X, y, sample_weight, seed=42)
        ),(
            CSRDataset32(X_csr32.data, X_csr32.indptr, X_csr32.indices, y32,
                         sample_weight32, seed=42),
            CSRDataset64(X_csr.data, X_csr.indptr, X_csr.indices, y,
                         sample_weight, seed=42)
        )
    ],
    ids=['ArrayDataset', 'CSRDataset']
)
def test_fused_types_consistency(dataset32, dataset64):
    NUMBER_OF_RUNS = 5
    for _ in range(NUMBER_OF_RUNS):
        # next sample
        xi32, yi32, _, _ = dataset32._next_py()
        xi64, yi64, _, _ = dataset64._next_py()

        xi_data32, _, _ = xi32
        xi_data64, _, _ = xi64

        assert xi_data32.dtype == np.float32
        assert xi_data64.dtype == np.float64

        assert_allclose(xi_data64, xi_data32, rtol=1e-5)
        assert_allclose(yi64, yi32, rtol=1e-5)
