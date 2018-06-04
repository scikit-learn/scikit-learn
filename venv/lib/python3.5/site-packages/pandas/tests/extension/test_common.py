import numpy as np
import pytest

import pandas as pd
import pandas.util.testing as tm
from pandas.core.arrays import ExtensionArray
from pandas.core.dtypes.common import is_extension_array_dtype
from pandas.core.dtypes import dtypes


class DummyDtype(dtypes.ExtensionDtype):
    pass


class DummyArray(ExtensionArray):

    def __init__(self, data):
        self.data = data

    def __array__(self, dtype):
        return self.data

    @property
    def dtype(self):
        return self.data.dtype


class TestExtensionArrayDtype(object):

    @pytest.mark.parametrize('values', [
        pd.Categorical([]),
        pd.Categorical([]).dtype,
        pd.Series(pd.Categorical([])),
        DummyDtype(),
        DummyArray(np.array([1, 2])),
    ])
    def test_is_extension_array_dtype(self, values):
        assert is_extension_array_dtype(values)

    @pytest.mark.parametrize('values', [
        np.array([]),
        pd.Series(np.array([])),
    ])
    def test_is_not_extension_array_dtype(self, values):
        assert not is_extension_array_dtype(values)


def test_astype():

    arr = DummyArray(np.array([1, 2, 3]))
    expected = np.array([1, 2, 3], dtype=object)

    result = arr.astype(object)
    tm.assert_numpy_array_equal(result, expected)

    result = arr.astype('object')
    tm.assert_numpy_array_equal(result, expected)


def test_astype_no_copy():
    arr = DummyArray(np.array([1, 2, 3], dtype=np.int64))
    result = arr.astype(arr.dtype, copy=False)

    assert arr.data is result

    result = arr.astype(arr.dtype)
    assert arr.data is not result


@pytest.mark.parametrize('dtype', [
    dtypes.DatetimeTZDtype('ns', 'US/Central'),
    dtypes.PeriodDtype("D"),
    dtypes.IntervalDtype(),
])
def test_is_not_extension_array_dtype(dtype):
    assert not isinstance(dtype, dtypes.ExtensionDtype)
    assert not is_extension_array_dtype(dtype)


@pytest.mark.parametrize('dtype', [
    dtypes.CategoricalDtype(),
])
def test_is_extension_array_dtype(dtype):
    assert isinstance(dtype, dtypes.ExtensionDtype)
    assert is_extension_array_dtype(dtype)
