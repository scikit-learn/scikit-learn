import numpy as np

import pandas as pd
from pandas.compat import StringIO
from pandas.core.dtypes.common import is_extension_array_dtype
from pandas.core.dtypes.dtypes import ExtensionDtype

from .base import BaseExtensionTests


class BaseInterfaceTests(BaseExtensionTests):
    """Tests that the basic interface is satisfied."""
    # ------------------------------------------------------------------------
    # Interface
    # ------------------------------------------------------------------------

    def test_len(self, data):
        assert len(data) == 100

    def test_ndim(self, data):
        assert data.ndim == 1

    def test_can_hold_na_valid(self, data):
        # GH-20761
        assert data._can_hold_na is True

    def test_memory_usage(self, data):
        s = pd.Series(data)
        result = s.memory_usage(index=False)
        assert result == s.nbytes

    def test_array_interface(self, data):
        result = np.array(data)
        assert result[0] == data[0]

    def test_repr(self, data):
        ser = pd.Series(data)
        assert data.dtype.name in repr(ser)

        df = pd.DataFrame({"A": data})
        repr(df)

    def test_dtype_name_in_info(self, data):
        buf = StringIO()
        pd.DataFrame({"A": data}).info(buf=buf)
        result = buf.getvalue()
        assert data.dtype.name in result

    def test_is_extension_array_dtype(self, data):
        assert is_extension_array_dtype(data)
        assert is_extension_array_dtype(data.dtype)
        assert is_extension_array_dtype(pd.Series(data))
        assert isinstance(data.dtype, ExtensionDtype)

    def test_no_values_attribute(self, data):
        # GH-20735: EA's with .values attribute give problems with internal
        # code, disallowing this for now until solved
        assert not hasattr(data, 'values')
        assert not hasattr(data, '_values')
