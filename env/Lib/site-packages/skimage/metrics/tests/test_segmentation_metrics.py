import numpy as np
import scipy as sp
import pytest

from skimage.metrics import (
    adapted_rand_error,
    variation_of_information,
    contingency_table,
)

from skimage._shared.testing import (
    assert_equal,
    assert_almost_equal,
    assert_array_equal,
)


@pytest.mark.parametrize("sparse_type", ["matrix", "array"])
def test_contingency_table(sparse_type):
    im_true = np.array([1, 2, 3, 4])
    im_test = np.array([1, 1, 8, 8])

    table1 = np.array(
        [
            [0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0],
            [0.0, 0.25, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0],
            [0.0, 0.25, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0],
            [0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.25],
            [0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.25],
        ]
    )

    sparse_table2 = contingency_table(
        im_true, im_test, normalize=True, sparse_type=sparse_type
    )
    table2 = sparse_table2.toarray()
    assert_array_equal(table1, table2)


def test_contingency_table_sparse_type():
    im_true = np.array([1, 2, 3, 4])
    im_test = np.array([1, 1, 8, 8])

    result = contingency_table(im_true, im_test)
    assert isinstance(result, sp.sparse.csr_matrix)

    result = contingency_table(im_true, im_test, sparse_type="matrix")
    assert isinstance(result, sp.sparse.csr_matrix)

    result = contingency_table(im_true, im_test, sparse_type="array")
    assert isinstance(result, sp.sparse.csr_array)

    with pytest.raises(ValueError, match="`sparse_type` must be 'array' or 'matrix'"):
        contingency_table(im_true, im_test, sparse_type="unknown")


def test_vi():
    im_true = np.array([1, 2, 3, 4])
    im_test = np.array([1, 1, 8, 8])
    assert_equal(np.sum(variation_of_information(im_true, im_test)), 1)


def test_vi_ignore_labels():
    im1 = np.array([[1, 0], [2, 3]], dtype='uint8')
    im2 = np.array([[1, 1], [1, 0]], dtype='uint8')

    false_splits, false_merges = variation_of_information(im1, im2, ignore_labels=[0])
    assert (false_splits, false_merges) == (0, 2 / 3)


def test_are():
    im_true = np.array([[2, 1], [1, 2]])
    im_test = np.array([[1, 2], [3, 1]])
    assert_almost_equal(adapted_rand_error(im_true, im_test), (0.3333333, 0.5, 1.0))
    assert_almost_equal(adapted_rand_error(im_true, im_test, alpha=0), (0, 0.5, 1.0))
    assert_almost_equal(adapted_rand_error(im_true, im_test, alpha=1), (0.5, 0.5, 1.0))

    with pytest.raises(ValueError):
        adapted_rand_error(im_true, im_test, alpha=1.01)
    with pytest.raises(ValueError):
        adapted_rand_error(im_true, im_test, alpha=-0.01)
