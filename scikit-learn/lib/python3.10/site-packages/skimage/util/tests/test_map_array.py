import numpy as np
import pytest
from skimage.util._map_array import map_array, ArrayMap

from skimage._shared import testing


_map_array_dtypes_in = [
    np.uint8,
    np.uint16,
    np.uint32,
    np.uint64,
    np.int8,
    np.int16,
    np.int32,
    np.int64,
]

_map_array_dtypes_out = _map_array_dtypes_in + [np.float32, np.float64]


@pytest.mark.parametrize("dtype_in", _map_array_dtypes_in)
@pytest.mark.parametrize("dtype_out", _map_array_dtypes_out)
@pytest.mark.parametrize("out_array", [True, False])
def test_map_array_simple(dtype_in, dtype_out, out_array):
    input_arr = np.array([0, 2, 0, 3, 4, 5, 0], dtype=dtype_in)
    input_vals = np.array([1, 2, 3, 4, 6], dtype=dtype_in)[::-1]
    output_vals = np.array([6, 7, 8, 9, 10], dtype=dtype_out)[::-1]
    desired = np.array([0, 7, 0, 8, 9, 0, 0], dtype=dtype_out)
    out = None
    if out_array:
        out = np.full(desired.shape, 11, dtype=dtype_out)
    result = map_array(
        input_arr=input_arr, input_vals=input_vals, output_vals=output_vals, out=out
    )
    np.testing.assert_array_equal(result, desired)
    assert result.dtype == dtype_out
    if out_array:
        assert out is result


def test_map_array_incorrect_output_shape():
    labels = np.random.randint(0, 5, size=(24, 25))
    out = np.empty((24, 24))
    in_values = np.unique(labels)
    out_values = np.random.random(in_values.shape).astype(out.dtype)
    with testing.raises(ValueError):
        map_array(labels, in_values, out_values, out=out)


def test_map_array_non_contiguous_output_array():
    labels = np.random.randint(0, 5, size=(24, 25))
    out = np.empty((24 * 3, 25 * 2))[::3, ::2]
    in_values = np.unique(labels)
    out_values = np.random.random(in_values.shape).astype(out.dtype)
    with testing.raises(ValueError):
        map_array(labels, in_values, out_values, out=out)


def test_arraymap_long_str():
    labels = np.random.randint(0, 40, size=(24, 25))
    in_values = np.unique(labels)
    out_values = np.random.random(in_values.shape)
    m = ArrayMap(in_values, out_values)
    assert len(str(m).split('\n')) == m._max_str_lines + 2


def test_arraymap_update():
    in_values = np.unique(np.random.randint(0, 200, size=5))
    out_values = np.random.random(len(in_values))
    m = ArrayMap(in_values, out_values)
    image = np.random.randint(1, len(m), size=(512, 512))
    assert np.all(m[image] < 1)  # missing values map to 0.
    m[1:] += 1
    assert np.all(m[image] >= 1)


def test_arraymap_bool_index():
    in_values = np.unique(np.random.randint(0, 200, size=5))
    out_values = np.random.random(len(in_values))
    m = ArrayMap(in_values, out_values)
    image = np.random.randint(1, len(in_values), size=(512, 512))
    assert np.all(m[image] < 1)  # missing values map to 0.
    positive = np.ones(len(m), dtype=bool)
    positive[0] = False
    m[positive] += 1
    assert np.all(m[image] >= 1)
