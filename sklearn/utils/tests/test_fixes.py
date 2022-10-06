# Authors: Gael Varoquaux <gael.varoquaux@normalesup.org>
#          Justin Vincent
#          Lars Buitinck
# License: BSD 3 clause

import math

import numpy as np
import pytest
import scipy.stats

from sklearn.utils._testing import assert_array_equal

from sklearn.utils.fixes import _object_dtype_isnan, loguniform, csr_hstack
from scipy.sparse import random as sparse_random, csr_matrix, hstack


@pytest.mark.parametrize("dtype, val", ([object, 1], [object, "a"], [float, 1]))
def test_object_dtype_isnan(dtype, val):
    X = np.array([[val, np.nan], [np.nan, val]], dtype=dtype)

    expected_mask = np.array([[False, True], [True, False]])

    mask = _object_dtype_isnan(X)

    assert_array_equal(mask, expected_mask)


@pytest.mark.parametrize("low,high,base", [(-1, 0, 10), (0, 2, np.exp(1)), (-1, 1, 2)])
def test_loguniform(low, high, base):
    rv = loguniform(base**low, base**high)
    assert isinstance(rv, scipy.stats._distn_infrastructure.rv_frozen)
    rvs = rv.rvs(size=2000, random_state=0)

    # Test the basics; right bounds, right size
    assert (base**low <= rvs).all() and (rvs <= base**high).all()
    assert len(rvs) == 2000

    # Test that it's actually (fairly) uniform
    log_rvs = np.array([math.log(x, base) for x in rvs])
    counts, _ = np.histogram(log_rvs)
    assert counts.mean() == 200
    assert np.abs(counts - counts.mean()).max() <= 40

    # Test that random_state works
    assert loguniform(base**low, base**high).rvs(random_state=0) == loguniform(
        base**low, base**high
    ).rvs(random_state=0)


# TODO: Remove when SciPy 1.10 is the minimum supported version
def test_csr_hstack():
    n_rows = 10
    msg = "No matrices were provided to stack"
    with pytest.raises(ValueError, match=msg):
        csr_hstack([])

    X1 = sparse_random(n_rows, 2, format="csr")
    X2 = sparse_random(n_rows + 1, 2, format="csr")
    msg = "Mismatching dimensions along axis*"
    with pytest.raises(ValueError, match=msg):
        csr_hstack([X1, X2])

    X1 = sparse_random(n_rows, 2, density=0, format="csr")
    X2 = sparse_random(n_rows, 2, density=0, format="csr")
    X_stacked = csr_hstack([X1, X2])
    assert X_stacked.data.size == 0
    assert X_stacked.indices.size == 0
    assert_array_equal(X_stacked.indptr, np.zeros(n_rows + 1, dtype=np.int64))


# TODO: Remove when SciPy 1.10 is the minimum supported version
def test_csr_hstack_int64():
    """
    Tests if hstack properly promotes to indices and indptr arrays to np.int64
    when using np.int32 during concatenation would result in either array
    overflowing.
    """
    max_int32 = np.iinfo(np.int32).max

    # First case: indices would overflow with int32
    data = [1.0]
    row = [0]

    max_indices_1 = max_int32 - 1
    max_indices_2 = 3

    # Individual indices arrays are representable with int32
    col_1 = [max_indices_1 - 1]
    col_2 = [max_indices_2 - 1]

    X_1 = csr_matrix((data, (row, col_1)))
    X_2 = csr_matrix((data, (row, col_2)))

    assert max(max_indices_1 - 1, max_indices_2 - 1) < max_int32
    assert X_1.indices.dtype == X_1.indptr.dtype == np.int32
    assert X_2.indices.dtype == X_2.indptr.dtype == np.int32

    # ... but when concatenating their CSR matrices, the resulting indices
    # array can't be represented with int32 and must be promoted to int64.
    X_hs = hstack([X_1, X_2], format="csr")

    assert X_hs.indices.max() == max_indices_1 + max_indices_2 - 1
    assert max_indices_1 + max_indices_2 - 1 > max_int32
    assert X_hs.indices.dtype == X_hs.indptr.dtype == np.int64

    # Even if the matrices are empty, we must account for their size
    # contribution so that we may safely set the final elements.
    X_1_empty = csr_matrix(X_1.shape)
    X_2_empty = csr_matrix(X_2.shape)
    X_hs_empty = hstack([X_1_empty, X_2_empty], format="csr")

    assert X_hs_empty.shape == X_hs.shape
    assert X_hs_empty.indices.dtype == np.int64

    # Should be just small enough to stay in int32 after stack. Note that
    # we theoretically could support indices.max() == max_int32, but due to an
    # edge-case in the underlying sparsetools code
    # (namely the `coo_tocsr` routine),
    # we require that max(X_hs_32.shape) < max_int32 as well.
    # Hence we can only support max_int32 - 1.
    col_3 = [max_int32 - max_indices_1 - 1]
    X_3 = csr_matrix((data, (row, col_3)))
    X_hs_32 = hstack([X_1, X_3], format="csr")
    assert X_hs_32.indices.dtype == np.int32
    assert X_hs_32.indices.max() == max_int32 - 1
