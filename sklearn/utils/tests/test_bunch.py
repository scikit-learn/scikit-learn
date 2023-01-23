import warnings

import numpy as np
import pytest

from sklearn.utils import Bunch


def test_bunch_attribute_deprecation():
    """Check that bunch raises deprecation message with `__getattr__`."""
    bunch = Bunch()
    values = np.asarray([1, 2, 3])
    msg = (
        "Key: 'values', is deprecated in 1.1 and will be "
        "removed in 1.3. Please use 'pdp_values' instead"
    )
    bunch._set_deprecated(
        values, new_key="pdp_values", deprecated_key="values", warning_message=msg
    )

    with warnings.catch_warnings():
        # Does not warn for "pdp_values"
        warnings.simplefilter("error")
        v = bunch["pdp_values"]

    assert v is values

    with pytest.warns(FutureWarning, match=msg):
        # Warns for "values"
        v = bunch["values"]

    assert v is values
