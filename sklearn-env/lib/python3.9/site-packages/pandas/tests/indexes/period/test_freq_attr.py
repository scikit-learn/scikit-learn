import pytest

from pandas import (
    offsets,
    period_range,
)
import pandas._testing as tm


class TestFreq:
    def test_freq_setter_deprecated(self):
        # GH#20678
        idx = period_range("2018Q1", periods=4, freq="Q")

        # no warning for getter
        with tm.assert_produces_warning(None):
            idx.freq

        # warning for setter
        with pytest.raises(AttributeError, match="can't set attribute"):
            idx.freq = offsets.Day()
