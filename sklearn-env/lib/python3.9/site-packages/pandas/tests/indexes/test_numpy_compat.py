import numpy as np
import pytest

from pandas import (
    CategoricalIndex,
    DatetimeIndex,
    Index,
    PeriodIndex,
    TimedeltaIndex,
    isna,
)
import pandas._testing as tm
from pandas.core.api import (
    Float64Index,
    NumericIndex,
)
from pandas.core.arrays import BooleanArray
from pandas.core.indexes.datetimelike import DatetimeIndexOpsMixin


@pytest.mark.parametrize(
    "func",
    [
        np.exp,
        np.exp2,
        np.expm1,
        np.log,
        np.log2,
        np.log10,
        np.log1p,
        np.sqrt,
        np.sin,
        np.cos,
        np.tan,
        np.arcsin,
        np.arccos,
        np.arctan,
        np.sinh,
        np.cosh,
        np.tanh,
        np.arcsinh,
        np.arccosh,
        np.arctanh,
        np.deg2rad,
        np.rad2deg,
    ],
    ids=lambda x: x.__name__,
)
def test_numpy_ufuncs_basic(index, func):
    # test ufuncs of numpy, see:
    # https://numpy.org/doc/stable/reference/ufuncs.html

    if isinstance(index, DatetimeIndexOpsMixin):
        with tm.external_error_raised((TypeError, AttributeError)):
            with np.errstate(all="ignore"):
                func(index)
    elif isinstance(index, NumericIndex) or (
        not isinstance(index.dtype, np.dtype) and index.dtype._is_numeric
    ):
        # coerces to float (e.g. np.sin)
        with np.errstate(all="ignore"):
            result = func(index)
            exp = Index(func(index.values), name=index.name)

        tm.assert_index_equal(result, exp)
        if type(index) is not Index:
            # i.e NumericIndex
            assert isinstance(result, Float64Index)
        else:
            # e.g. np.exp with Int64 -> Float64
            assert type(result) is Index
    else:
        # raise AttributeError or TypeError
        if len(index) == 0:
            pass
        else:
            with tm.external_error_raised((TypeError, AttributeError)):
                with np.errstate(all="ignore"):
                    func(index)


@pytest.mark.parametrize(
    "func", [np.isfinite, np.isinf, np.isnan, np.signbit], ids=lambda x: x.__name__
)
def test_numpy_ufuncs_other(index, func, request):
    # test ufuncs of numpy, see:
    # https://numpy.org/doc/stable/reference/ufuncs.html
    if isinstance(index, (DatetimeIndex, TimedeltaIndex)):

        if func in (np.isfinite, np.isinf, np.isnan):
            # numpy 1.18 changed isinf and isnan to not raise on dt64/td64
            result = func(index)
            assert isinstance(result, np.ndarray)
        else:
            with tm.external_error_raised(TypeError):
                func(index)

    elif isinstance(index, PeriodIndex):
        with tm.external_error_raised(TypeError):
            func(index)

    elif isinstance(index, NumericIndex) or (
        not isinstance(index.dtype, np.dtype) and index.dtype._is_numeric
    ):
        # Results in bool array
        result = func(index)
        if not isinstance(index.dtype, np.dtype):
            # e.g. Int64 we expect to get BooleanArray back
            assert isinstance(result, BooleanArray)
        else:
            assert isinstance(result, np.ndarray)
        assert not isinstance(result, Index)
    else:
        if len(index) == 0:
            pass
        else:
            with tm.external_error_raised(TypeError):
                func(index)


@pytest.mark.parametrize("func", [np.maximum, np.minimum])
def test_numpy_ufuncs_reductions(index, func, request):
    # TODO: overlap with tests.series.test_ufunc.test_reductions
    if len(index) == 0:
        return

    if repr(index.dtype) == "string[pyarrow]":
        mark = pytest.mark.xfail(reason="ArrowStringArray has no min/max")
        request.node.add_marker(mark)

    if isinstance(index, CategoricalIndex) and index.dtype.ordered is False:
        with pytest.raises(TypeError, match="is not ordered for"):
            func.reduce(index)
        return
    else:
        result = func.reduce(index)

    if func is np.maximum:
        expected = index.max(skipna=False)
    else:
        expected = index.min(skipna=False)
        # TODO: do we have cases both with and without NAs?

    assert type(result) is type(expected)
    if isna(result):
        assert isna(expected)
    else:
        assert result == expected
