import numpy as np

import pytest

from sklearn.utils._readonly_array_wrapper import ReadonlyWrapper, _test_sum
from sklearn.utils._testing import create_memmap_backed_data


@pytest.mark.parametrize("readonly", ["flag", "memmap"])
@pytest.mark.parametrize("dtype", [np.float32, np.float64, np.int32, np.int64])
def test_readonly_array_wrapper(readonly, dtype):
    """Test that ReadonlyWrapper works as expected."""
    x = np.arange(10).astype(dtype)
    sum_origin = _test_sum(x)

    if readonly == "flag":
        x_readonly = x.copy()
        x_readonly.flags["WRITEABLE"] = False
    else:
        x_readonly = create_memmap_backed_data(x)

    with pytest.raises(ValueError, match="buffer source array is read-only"):
        _test_sum(x_readonly)

    x_readonly = ReadonlyWrapper(x_readonly)
    sum_readonly = _test_sum(x_readonly)
    assert sum_readonly == pytest.approx(sum_origin, rel=1e-11)
