import pytest
import numpy as np

import dask.array as da
from dask.array.numpy_compat import _make_sliced_dtype
from dask.array.utils import assert_eq


@pytest.fixture(params=[
    [('A', ('f4', (3, 2))), ('B', ('f4', 3)), ('C', ('f8', 3))],
    [('A', ('i4', (3, 2))), ('B', ('f4', 3)), ('C', ('S4', 3))],
])
def dtype(request):
    return np.dtype(request.param)


@pytest.fixture(params=[
    ['A'],
    ['A', 'B'],
    ['A', 'B', 'C'],
])
def index(request):
    return request.param


def test_basic():
    # sanity check
    dtype = [('a', 'f8'), ('b', 'f8'), ('c', 'f8')]
    x = np.ones((5, 3), dtype=dtype)
    dx = da.ones((5, 3), dtype=dtype, chunks=3)
    result = dx[['a', 'b']]
    expected = x[['a', 'b']]
    assert_eq(result, expected)


def test_slice_dtype(dtype, index):
    result = _make_sliced_dtype(dtype, index)
    expected = np.ones((5, len(dtype)), dtype=dtype)[index].dtype
    assert result == expected
