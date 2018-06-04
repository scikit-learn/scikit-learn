import pandas as pd
from pandas.core.internals import ObjectBlock

from .base import BaseExtensionTests


class BaseCastingTests(BaseExtensionTests):
    """Casting to and from ExtensionDtypes"""

    def test_astype_object_series(self, all_data):
        ser = pd.Series({"A": all_data})
        result = ser.astype(object)
        assert isinstance(result._data.blocks[0], ObjectBlock)

    def test_tolist(self, data):
        result = pd.Series(data).tolist()
        expected = list(data)
        assert result == expected

    def test_astype_str(self, data):
        result = pd.Series(data[:5]).astype(str)
        expected = pd.Series(data[:5].astype(str))
        self.assert_series_equal(result, expected)
