import numpy as np
import pytest

from sklearn.utils._protocols import DataFrameInterchangeProtocol
from sklearn.utils._testing import _convert_container


@pytest.mark.parametrize(
    "constructor_name, minversion",
    [
        ("dataframe", "1.5.0"),
        ("pyarrow", "12.0.0"),
        ("polars", "0.18.2"),
    ],
)
def test_dataframe_interchange_protocol(constructor_name, minversion):
    """Check that the protocol works with isinstance."""
    X = np.asarray([[1, 2, 3], [3, 4, 5]])
    columns_name = ["a", "b", "c"]
    df = _convert_container(X, constructor_name, columns_name, minversion=minversion)

    assert isinstance(df.__dataframe__(), DataFrameInterchangeProtocol)
