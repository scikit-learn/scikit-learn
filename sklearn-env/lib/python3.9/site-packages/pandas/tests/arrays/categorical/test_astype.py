import numpy as np
import pytest

from pandas import (
    Categorical,
    CategoricalDtype,
    NaT,
    Timestamp,
    array,
    to_datetime,
)
import pandas._testing as tm


class TestAstype:
    def test_astype_str_int_categories_to_nullable_int(self):
        # GH#39616
        dtype = CategoricalDtype([str(i) for i in range(5)])
        codes = np.random.randint(5, size=20)
        arr = Categorical.from_codes(codes, dtype=dtype)

        res = arr.astype("Int64")
        expected = array(codes, dtype="Int64")
        tm.assert_extension_array_equal(res, expected)

    @pytest.mark.parametrize("ordered", [True, False])
    def test_astype(self, ordered):
        # string
        cat = Categorical(list("abbaaccc"), ordered=ordered)
        result = cat.astype(object)
        expected = np.array(cat)
        tm.assert_numpy_array_equal(result, expected)

        msg = r"Cannot cast object dtype to float64"
        with pytest.raises(ValueError, match=msg):
            cat.astype(float)

        # numeric
        cat = Categorical([0, 1, 2, 2, 1, 0, 1, 0, 2], ordered=ordered)
        result = cat.astype(object)
        expected = np.array(cat, dtype=object)
        tm.assert_numpy_array_equal(result, expected)

        result = cat.astype(int)
        expected = np.array(cat, dtype="int")
        tm.assert_numpy_array_equal(result, expected)

        result = cat.astype(float)
        expected = np.array(cat, dtype=float)
        tm.assert_numpy_array_equal(result, expected)

    @pytest.mark.parametrize("dtype_ordered", [True, False])
    @pytest.mark.parametrize("cat_ordered", [True, False])
    def test_astype_category(self, dtype_ordered, cat_ordered):
        # GH#10696/GH#18593
        data = list("abcaacbab")
        cat = Categorical(data, categories=list("bac"), ordered=cat_ordered)

        # standard categories
        dtype = CategoricalDtype(ordered=dtype_ordered)
        result = cat.astype(dtype)
        expected = Categorical(data, categories=cat.categories, ordered=dtype_ordered)
        tm.assert_categorical_equal(result, expected)

        # non-standard categories
        dtype = CategoricalDtype(list("adc"), dtype_ordered)
        result = cat.astype(dtype)
        expected = Categorical(data, dtype=dtype)
        tm.assert_categorical_equal(result, expected)

        if dtype_ordered is False:
            # dtype='category' can't specify ordered, so only test once
            result = cat.astype("category")
            expected = cat
            tm.assert_categorical_equal(result, expected)

    def test_astype_object_datetime_categories(self):
        # GH#40754
        cat = Categorical(to_datetime(["2021-03-27", NaT]))
        result = cat.astype(object)
        expected = np.array([Timestamp("2021-03-27 00:00:00"), NaT], dtype="object")
        tm.assert_numpy_array_equal(result, expected)

    def test_astype_object_timestamp_categories(self):
        # GH#18024
        cat = Categorical([Timestamp("2014-01-01")])
        result = cat.astype(object)
        expected = np.array([Timestamp("2014-01-01 00:00:00")], dtype="object")
        tm.assert_numpy_array_equal(result, expected)
