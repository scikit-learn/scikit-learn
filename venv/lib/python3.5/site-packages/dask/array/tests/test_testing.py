import sys

import pytest
import numpy as np

import dask.array as da
from dask.array.utils import assert_eq


@pytest.mark.skipif(sys.flags.optimize, reason="Assertions disabled.")
def test_assert_eq_checks_scalars():
    # https://github.com/dask/dask/issues/2680
    with pytest.raises(AssertionError):
        assert_eq(np.array(0), np.array(1))

    a = da.from_array(np.array([0]), 1)[0]
    b = np.array([1])[0]
    with pytest.raises(AssertionError):
        assert_eq(a, b)
