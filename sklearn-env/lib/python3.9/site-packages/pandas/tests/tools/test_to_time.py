from datetime import time
import locale

import numpy as np
import pytest

from pandas import Series
import pandas._testing as tm
from pandas.core.tools.datetimes import to_time as to_time_alias
from pandas.core.tools.times import to_time

fails_on_non_english = pytest.mark.xfail(
    locale.getlocale()[0] in ("zh_CN", "it_IT"),
    reason="fail on a CI build with LC_ALL=zh_CN.utf8/it_IT.utf8",
)


class TestToTime:
    @pytest.mark.parametrize(
        "time_string",
        [
            "14:15",
            "1415",
            pytest.param("2:15pm", marks=fails_on_non_english),
            pytest.param("0215pm", marks=fails_on_non_english),
            "14:15:00",
            "141500",
            pytest.param("2:15:00pm", marks=fails_on_non_english),
            pytest.param("021500pm", marks=fails_on_non_english),
            time(14, 15),
        ],
    )
    def test_parsers_time(self, time_string):
        # GH#11818
        assert to_time(time_string) == time(14, 15)

    def test_odd_format(self):
        new_string = "14.15"
        msg = r"Cannot convert arg \['14\.15'\] to a time"
        with pytest.raises(ValueError, match=msg):
            to_time(new_string)
        assert to_time(new_string, format="%H.%M") == time(14, 15)

    def test_arraylike(self):
        arg = ["14:15", "20:20"]
        expected_arr = [time(14, 15), time(20, 20)]
        assert to_time(arg) == expected_arr
        assert to_time(arg, format="%H:%M") == expected_arr
        assert to_time(arg, infer_time_format=True) == expected_arr
        assert to_time(arg, format="%I:%M%p", errors="coerce") == [None, None]

        res = to_time(arg, format="%I:%M%p", errors="ignore")
        tm.assert_numpy_array_equal(res, np.array(arg, dtype=np.object_))

        msg = "Cannot convert.+to a time with given format"
        with pytest.raises(ValueError, match=msg):
            to_time(arg, format="%I:%M%p", errors="raise")

        tm.assert_series_equal(
            to_time(Series(arg, name="test")), Series(expected_arr, name="test")
        )

        res = to_time(np.array(arg))
        assert isinstance(res, list)
        assert res == expected_arr


def test_to_time_alias():
    expected = time(14, 15)

    with tm.assert_produces_warning(FutureWarning):
        result = to_time_alias(expected)

    assert result == expected
