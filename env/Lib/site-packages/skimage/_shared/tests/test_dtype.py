import numpy as np
import pytest


from ..dtype import numeric_dtype_min_max, numeric_types


class Test_numeric_dtype_min_max:
    @pytest.mark.parametrize("dtype", numeric_types)
    def test_all_numeric_types(self, dtype):
        min_, max_ = numeric_dtype_min_max(dtype)
        assert np.isscalar(min_)
        assert np.isscalar(max_)
        assert min_ < max_
