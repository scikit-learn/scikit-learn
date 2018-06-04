import decimal

import numpy as np
import pandas as pd
import pandas.util.testing as tm
import pytest

from pandas.tests.extension import base

from .array import DecimalDtype, DecimalArray, make_data


@pytest.fixture
def dtype():
    return DecimalDtype()


@pytest.fixture
def data():
    return DecimalArray(make_data())


@pytest.fixture
def data_missing():
    return DecimalArray([decimal.Decimal('NaN'), decimal.Decimal(1)])


@pytest.fixture
def data_for_sorting():
    return DecimalArray([decimal.Decimal('1'),
                         decimal.Decimal('2'),
                         decimal.Decimal('0')])


@pytest.fixture
def data_missing_for_sorting():
    return DecimalArray([decimal.Decimal('1'),
                         decimal.Decimal('NaN'),
                         decimal.Decimal('0')])


@pytest.fixture
def na_cmp():
    return lambda x, y: x.is_nan() and y.is_nan()


@pytest.fixture
def na_value():
    return decimal.Decimal("NaN")


@pytest.fixture
def data_for_grouping():
    b = decimal.Decimal('1.0')
    a = decimal.Decimal('0.0')
    c = decimal.Decimal('2.0')
    na = decimal.Decimal('NaN')
    return DecimalArray([b, b, na, na, a, a, b, c])


class BaseDecimal(object):

    def assert_series_equal(self, left, right, *args, **kwargs):

        left_na = left.isna()
        right_na = right.isna()

        tm.assert_series_equal(left_na, right_na)
        return tm.assert_series_equal(left[~left_na],
                                      right[~right_na],
                                      *args, **kwargs)

    def assert_frame_equal(self, left, right, *args, **kwargs):
        # TODO(EA): select_dtypes
        tm.assert_index_equal(
            left.columns, right.columns,
            exact=kwargs.get('check_column_type', 'equiv'),
            check_names=kwargs.get('check_names', True),
            check_exact=kwargs.get('check_exact', False),
            check_categorical=kwargs.get('check_categorical', True),
            obj='{obj}.columns'.format(obj=kwargs.get('obj', 'DataFrame')))

        decimals = (left.dtypes == 'decimal').index

        for col in decimals:
            self.assert_series_equal(left[col], right[col],
                                     *args, **kwargs)

        left = left.drop(columns=decimals)
        right = right.drop(columns=decimals)
        tm.assert_frame_equal(left, right, *args, **kwargs)


class TestDtype(BaseDecimal, base.BaseDtypeTests):
    pass


class TestInterface(BaseDecimal, base.BaseInterfaceTests):
    pass


class TestConstructors(BaseDecimal, base.BaseConstructorsTests):
    pass


class TestReshaping(BaseDecimal, base.BaseReshapingTests):
    pass


class TestGetitem(BaseDecimal, base.BaseGetitemTests):

    def test_take_na_value_other_decimal(self):
        arr = DecimalArray([decimal.Decimal('1.0'),
                            decimal.Decimal('2.0')])
        result = arr.take([0, -1], allow_fill=True,
                          fill_value=decimal.Decimal('-1.0'))
        expected = DecimalArray([decimal.Decimal('1.0'),
                                 decimal.Decimal('-1.0')])
        self.assert_extension_array_equal(result, expected)


class TestMissing(BaseDecimal, base.BaseMissingTests):
    pass


class TestMethods(BaseDecimal, base.BaseMethodsTests):
    @pytest.mark.parametrize('dropna', [True, False])
    @pytest.mark.xfail(reason="value_counts not implemented yet.")
    def test_value_counts(self, all_data, dropna):
        all_data = all_data[:10]
        if dropna:
            other = np.array(all_data[~all_data.isna()])
        else:
            other = all_data

        result = pd.Series(all_data).value_counts(dropna=dropna).sort_index()
        expected = pd.Series(other).value_counts(dropna=dropna).sort_index()

        tm.assert_series_equal(result, expected)


class TestCasting(BaseDecimal, base.BaseCastingTests):
    pass


class TestGroupby(BaseDecimal, base.BaseGroupbyTests):
    pass


def test_series_constructor_coerce_data_to_extension_dtype_raises():
    xpr = ("Cannot cast data to extension dtype 'decimal'. Pass the "
           "extension array directly.")
    with tm.assert_raises_regex(ValueError, xpr):
        pd.Series([0, 1, 2], dtype=DecimalDtype())


def test_series_constructor_with_same_dtype_ok():
    arr = DecimalArray([decimal.Decimal('10.0')])
    result = pd.Series(arr, dtype=DecimalDtype())
    expected = pd.Series(arr)
    tm.assert_series_equal(result, expected)


def test_series_constructor_coerce_extension_array_to_dtype_raises():
    arr = DecimalArray([decimal.Decimal('10.0')])
    xpr = r"Cannot specify a dtype 'int64' .* \('decimal'\)."

    with tm.assert_raises_regex(ValueError, xpr):
        pd.Series(arr, dtype='int64')


def test_dataframe_constructor_with_same_dtype_ok():
    arr = DecimalArray([decimal.Decimal('10.0')])

    result = pd.DataFrame({"A": arr}, dtype=DecimalDtype())
    expected = pd.DataFrame({"A": arr})
    tm.assert_frame_equal(result, expected)


def test_dataframe_constructor_with_different_dtype_raises():
    arr = DecimalArray([decimal.Decimal('10.0')])

    xpr = "Cannot coerce extension array to dtype 'int64'. "
    with tm.assert_raises_regex(ValueError, xpr):
        pd.DataFrame({"A": arr}, dtype='int64')
