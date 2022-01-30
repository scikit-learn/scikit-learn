from __future__ import annotations

import numpy as np
import pytest

import pandas as pd
import pandas._testing as tm
from pandas.core import ops
from pandas.tests.extension.base.base import BaseExtensionTests


class BaseOpsUtil(BaseExtensionTests):
    def get_op_from_name(self, op_name: str):
        return tm.get_op_from_name(op_name)

    def check_opname(self, ser: pd.Series, op_name: str, other, exc=Exception):
        op = self.get_op_from_name(op_name)

        self._check_op(ser, op, other, op_name, exc)

    def _combine(self, obj, other, op):
        if isinstance(obj, pd.DataFrame):
            if len(obj.columns) != 1:
                raise NotImplementedError
            expected = obj.iloc[:, 0].combine(other, op).to_frame()
        else:
            expected = obj.combine(other, op)
        return expected

    def _check_op(
        self, ser: pd.Series, op, other, op_name: str, exc=NotImplementedError
    ):
        if exc is None:
            result = op(ser, other)
            expected = self._combine(ser, other, op)
            assert isinstance(result, type(ser))
            self.assert_equal(result, expected)
        else:
            with pytest.raises(exc):
                op(ser, other)

    def _check_divmod_op(self, ser: pd.Series, op, other, exc=Exception):
        # divmod has multiple return values, so check separately
        if exc is None:
            result_div, result_mod = op(ser, other)
            if op is divmod:
                expected_div, expected_mod = ser // other, ser % other
            else:
                expected_div, expected_mod = other // ser, other % ser
            self.assert_series_equal(result_div, expected_div)
            self.assert_series_equal(result_mod, expected_mod)
        else:
            with pytest.raises(exc):
                divmod(ser, other)


class BaseArithmeticOpsTests(BaseOpsUtil):
    """
    Various Series and DataFrame arithmetic ops methods.

    Subclasses supporting various ops should set the class variables
    to indicate that they support ops of that kind

    * series_scalar_exc = TypeError
    * frame_scalar_exc = TypeError
    * series_array_exc = TypeError
    * divmod_exc = TypeError
    """

    series_scalar_exc: type[TypeError] | None = TypeError
    frame_scalar_exc: type[TypeError] | None = TypeError
    series_array_exc: type[TypeError] | None = TypeError
    divmod_exc: type[TypeError] | None = TypeError

    def test_arith_series_with_scalar(self, data, all_arithmetic_operators):
        # series & scalar
        op_name = all_arithmetic_operators
        ser = pd.Series(data)
        self.check_opname(ser, op_name, ser.iloc[0], exc=self.series_scalar_exc)

    def test_arith_frame_with_scalar(self, data, all_arithmetic_operators):
        # frame & scalar
        op_name = all_arithmetic_operators
        df = pd.DataFrame({"A": data})
        self.check_opname(df, op_name, data[0], exc=self.frame_scalar_exc)

    def test_arith_series_with_array(self, data, all_arithmetic_operators):
        # ndarray & other series
        op_name = all_arithmetic_operators
        ser = pd.Series(data)
        self.check_opname(
            ser, op_name, pd.Series([ser.iloc[0]] * len(ser)), exc=self.series_array_exc
        )

    def test_divmod(self, data):
        ser = pd.Series(data)
        self._check_divmod_op(ser, divmod, 1, exc=self.divmod_exc)
        self._check_divmod_op(1, ops.rdivmod, ser, exc=self.divmod_exc)

    def test_divmod_series_array(self, data, data_for_twos):
        ser = pd.Series(data)
        self._check_divmod_op(ser, divmod, data)

        other = data_for_twos
        self._check_divmod_op(other, ops.rdivmod, ser)

        other = pd.Series(other)
        self._check_divmod_op(other, ops.rdivmod, ser)

    def test_add_series_with_extension_array(self, data):
        ser = pd.Series(data)
        result = ser + data
        expected = pd.Series(data + data)
        self.assert_series_equal(result, expected)

    @pytest.mark.parametrize("box", [pd.Series, pd.DataFrame])
    def test_direct_arith_with_ndframe_returns_not_implemented(self, data, box):
        # EAs should return NotImplemented for ops with Series/DataFrame
        # Pandas takes care of unboxing the series and calling the EA's op.
        other = pd.Series(data)
        if box is pd.DataFrame:
            other = other.to_frame()
        if hasattr(data, "__add__"):
            result = data.__add__(other)
            assert result is NotImplemented
        else:
            raise pytest.skip(f"{type(data).__name__} does not implement add")


class BaseComparisonOpsTests(BaseOpsUtil):
    """Various Series and DataFrame comparison ops methods."""

    def _compare_other(self, ser: pd.Series, data, op, other):

        if op.__name__ in ["eq", "ne"]:
            # comparison should match point-wise comparisons
            result = op(ser, other)
            expected = ser.combine(other, op)
            self.assert_series_equal(result, expected)

        else:
            exc = None
            try:
                result = op(ser, other)
            except Exception as err:
                exc = err

            if exc is None:
                # Didn't error, then should match pointwise behavior
                expected = ser.combine(other, op)
                self.assert_series_equal(result, expected)
            else:
                with pytest.raises(type(exc)):
                    ser.combine(other, op)

    def test_compare_scalar(self, data, comparison_op):
        ser = pd.Series(data)
        self._compare_other(ser, data, comparison_op, 0)

    def test_compare_array(self, data, comparison_op):
        ser = pd.Series(data)
        other = pd.Series([data[0]] * len(data))
        self._compare_other(ser, data, comparison_op, other)

    @pytest.mark.parametrize("box", [pd.Series, pd.DataFrame])
    def test_direct_arith_with_ndframe_returns_not_implemented(self, data, box):
        # EAs should return NotImplemented for ops with Series/DataFrame
        # Pandas takes care of unboxing the series and calling the EA's op.
        other = pd.Series(data)
        if box is pd.DataFrame:
            other = other.to_frame()

        if hasattr(data, "__eq__"):
            result = data.__eq__(other)
            assert result is NotImplemented
        else:
            raise pytest.skip(f"{type(data).__name__} does not implement __eq__")

        if hasattr(data, "__ne__"):
            result = data.__ne__(other)
            assert result is NotImplemented
        else:
            raise pytest.skip(f"{type(data).__name__} does not implement __ne__")


class BaseUnaryOpsTests(BaseOpsUtil):
    def test_invert(self, data):
        ser = pd.Series(data, name="name")
        result = ~ser
        expected = pd.Series(~data, name="name")
        self.assert_series_equal(result, expected)

    @pytest.mark.parametrize("ufunc", [np.positive, np.negative, np.abs])
    def test_unary_ufunc_dunder_equivalence(self, data, ufunc):
        # the dunder __pos__ works if and only if np.positive works,
        #  same for __neg__/np.negative and __abs__/np.abs
        attr = {np.positive: "__pos__", np.negative: "__neg__", np.abs: "__abs__"}[
            ufunc
        ]

        exc = None
        try:
            result = getattr(data, attr)()
        except Exception as err:
            exc = err

            # if __pos__ raised, then so should the ufunc
            with pytest.raises((type(exc), TypeError)):
                ufunc(data)
        else:
            alt = ufunc(data)
            self.assert_extension_array_equal(result, alt)
