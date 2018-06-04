from pandas.core.series import Series

import pytest
import pandas.util.testing as tm


@pytest.fixture
def series():
    return Series([1, 2, 3, 4, 5])


class TestSeriesValidate(object):
    """Tests for error handling related to data types of method arguments."""

    @pytest.mark.parametrize("func", ["reset_index", "_set_name",
                                      "sort_values", "sort_index",
                                      "rename", "dropna"])
    @pytest.mark.parametrize("inplace", [1, "True", [1, 2, 3], 5.0])
    def test_validate_bool_args(self, series, func, inplace):
        msg = "For argument \"inplace\" expected type bool"
        kwargs = dict(inplace=inplace)

        if func == "_set_name":
            kwargs["name"] = "hello"

        with tm.assert_raises_regex(ValueError, msg):
            getattr(series, func)(**kwargs)
