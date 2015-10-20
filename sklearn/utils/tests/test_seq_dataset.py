# Author: Tom Dupre la Tour <tom.dupre-la-tour@m4x.org>
#
# License: BSD 3 clause

import numpy as np
import scipy.sparse as sp

from sklearn.utils.seq_dataset import ArrayDataset, CSRDataset
from sklearn.datasets import load_iris

from numpy.testing import assert_array_equal
from nose.tools import assert_equal

iris = load_iris()
X = iris.data.astype(np.float64)
y = iris.target.astype(np.float64)
X_csr = sp.csr_matrix(X)
sample_weight = np.arange(y.size, dtype=np.float64)


def test_seq_dataset():
    dataset1 = ArrayDataset(X, y, sample_weight, seed=42)
    dataset2 = CSRDataset(X_csr.data, X_csr.indptr, X_csr.indices,
                          y, sample_weight, seed=42)

    for dataset in (dataset1, dataset2):
        for i in range(5):
            # next sample
            xi_, yi, swi, idx = dataset._next_py()
            xi = sp.csr_matrix((xi_), shape=(1, X.shape[1]))

            assert_array_equal(xi.data, X_csr[idx].data)
            assert_array_equal(xi.indices, X_csr[idx].indices)
            assert_array_equal(xi.indptr, X_csr[idx].indptr)
            assert_equal(yi, y[idx])
            assert_equal(swi, sample_weight[idx])

            # random sample
            xi_, yi, swi, idx = dataset._random_py()
            xi = sp.csr_matrix((xi_), shape=(1, X.shape[1]))

            assert_array_equal(xi.data, X_csr[idx].data)
            assert_array_equal(xi.indices, X_csr[idx].indices)
            assert_array_equal(xi.indptr, X_csr[idx].indptr)
            assert_equal(yi, y[idx])
            assert_equal(swi, sample_weight[idx])


def test_seq_dataset_shuffle():
    dataset1 = ArrayDataset(X, y, sample_weight, seed=42)
    dataset2 = CSRDataset(X_csr.data, X_csr.indptr, X_csr.indices,
                          y, sample_weight, seed=42)

    # not shuffled
    for i in range(5):
        _, _, _, idx1 = dataset1._next_py()
        _, _, _, idx2 = dataset2._next_py()
        assert_equal(idx1, i)
        assert_equal(idx2, i)

    for i in range(5):
        _, _, _, idx1 = dataset1._random_py()
        _, _, _, idx2 = dataset2._random_py()
        assert_equal(idx1, idx2)

    seed = 77
    dataset1._shuffle_py(seed)
    dataset2._shuffle_py(seed)

    for i in range(5):
        _, _, _, idx1 = dataset1._next_py()
        _, _, _, idx2 = dataset2._next_py()
        assert_equal(idx1, idx2)

        _, _, _, idx1 = dataset1._random_py()
        _, _, _, idx2 = dataset2._random_py()
        assert_equal(idx1, idx2)

