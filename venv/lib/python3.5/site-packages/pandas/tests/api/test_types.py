# -*- coding: utf-8 -*-

import pytest

from warnings import catch_warnings

import pandas
from pandas.api import types
from pandas.util import testing as tm

from .test_api import Base


class TestTypes(Base):

    allowed = ['is_bool', 'is_bool_dtype',
               'is_categorical', 'is_categorical_dtype', 'is_complex',
               'is_complex_dtype', 'is_datetime64_any_dtype',
               'is_datetime64_dtype', 'is_datetime64_ns_dtype',
               'is_datetime64tz_dtype', 'is_datetimetz', 'is_dtype_equal',
               'is_extension_type', 'is_float', 'is_float_dtype',
               'is_int64_dtype', 'is_integer',
               'is_integer_dtype', 'is_number', 'is_numeric_dtype',
               'is_object_dtype', 'is_scalar', 'is_sparse',
               'is_string_dtype', 'is_signed_integer_dtype',
               'is_timedelta64_dtype', 'is_timedelta64_ns_dtype',
               'is_unsigned_integer_dtype', 'is_period',
               'is_period_dtype', 'is_interval', 'is_interval_dtype',
               'is_re', 'is_re_compilable',
               'is_dict_like', 'is_iterator', 'is_file_like',
               'is_list_like', 'is_hashable', 'is_array_like',
               'is_named_tuple',
               'pandas_dtype', 'union_categoricals', 'infer_dtype']
    deprecated = ['is_any_int_dtype', 'is_floating_dtype', 'is_sequence']
    dtypes = ['CategoricalDtype', 'DatetimeTZDtype',
              'PeriodDtype', 'IntervalDtype']

    def test_types(self):

        self.check(types, self.allowed + self.dtypes + self.deprecated)

    def check_deprecation(self, fold, fnew):
        with tm.assert_produces_warning(DeprecationWarning):
            try:
                result = fold('foo')
                expected = fnew('foo')
                assert result == expected
            except TypeError:
                pytest.raises(TypeError, lambda: fnew('foo'))
            except AttributeError:
                pytest.raises(AttributeError, lambda: fnew('foo'))

    def test_deprecated_from_api_types(self):

        for t in self.deprecated:
            with tm.assert_produces_warning(FutureWarning,
                                            check_stacklevel=False):
                getattr(types, t)(1)


def test_moved_infer_dtype():

    with catch_warnings(record=True):
        e = pandas.lib.infer_dtype('foo')
        assert e is not None
