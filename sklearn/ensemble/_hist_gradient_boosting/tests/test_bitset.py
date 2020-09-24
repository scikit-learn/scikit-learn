import pytest
import numpy as np
from numpy.testing import assert_allclose

from sklearn.ensemble._hist_gradient_boosting._bitset import (
    set_bitset_memoryview,
    set_raw_bitset_from_binned_bitset
)


@pytest.mark.parametrize("expected_bitset, values_to_insert", [
    (np.array([2**0 + 2**4, 2**1, 0], dtype=np.uint32), [0, 4, 33]),
    (np.array([2**31, 2**0, 2**15], dtype=np.uint32), [31, 32, 79])
])
def test_set_bitset_memoryview(expected_bitset, values_to_insert):
    bitset = np.zeros(3, dtype=np.uint32)
    for value in values_to_insert:
        set_bitset_memoryview(bitset, value)
    assert_allclose(expected_bitset, bitset)


@pytest.mark.parametrize(
    "expected_raw_bitset, binned_bitset, categories", [
        (np.array([2**3 + 2**5 + 2**31, 2**0 + 2**11], dtype=np.uint32),
         np.array([2**0 + 2**2 + 2**4 + 2**5 + 2**6], dtype=np.uint32),
         np.array([3.0, 4.0, 5.0, 10.0, 31.0, 32.0, 43.0], dtype=np.float64)),
        (np.array([0, 2**1 + 2**20], dtype=np.uint32),
         np.array([2**1 + 2**3], dtype=np.uint32),
         np.array([3.0, 33.0, 50.0, 52.0], dtype=np.float64)),
        (np.array([2**20, 2**30], dtype=np.uint32),
         np.array([2**10 + 2**31], dtype=np.uint32),
         np.array(list(range(0, 64, 2)), dtype=np.float64))
    ])
def test_raw_bitset_mv(expected_raw_bitset, binned_bitset, categories):
    raw_bitset = np.zeros(2, dtype=np.uint32)
    set_raw_bitset_from_binned_bitset(raw_bitset, binned_bitset, categories)
    assert_allclose(expected_raw_bitset, raw_bitset)
