import numpy as np
import pandas as pd

from .base import BaseExtensionTests


class BaseDtypeTests(BaseExtensionTests):
    """Base class for ExtensionDtype classes"""

    def test_name(self, dtype):
        assert isinstance(dtype.name, str)

    def test_kind(self, dtype):
        valid = set('biufcmMOSUV')
        if dtype.kind is not None:
            assert dtype.kind in valid

    def test_construct_from_string_own_name(self, dtype):
        result = dtype.construct_from_string(dtype.name)
        assert type(result) is type(dtype)

        # check OK as classmethod
        result = type(dtype).construct_from_string(dtype.name)
        assert type(result) is type(dtype)

    def test_is_dtype_from_name(self, dtype):
        result = type(dtype).is_dtype(dtype.name)
        assert result is True

    def test_is_dtype_unboxes_dtype(self, data, dtype):
        assert dtype.is_dtype(data) is True

    def test_is_dtype_from_self(self, dtype):
        result = type(dtype).is_dtype(dtype)
        assert result is True

    def test_is_not_string_type(self, dtype):
        return not pd.api.types.is_string_dtype(dtype)

    def test_is_not_object_type(self, dtype):
        return not pd.api.types.is_object_dtype(dtype)

    def test_eq_with_str(self, dtype):
        assert dtype == dtype.name
        assert dtype != dtype.name + '-suffix'

    def test_eq_with_numpy_object(self, dtype):
        assert dtype != np.dtype('object')
