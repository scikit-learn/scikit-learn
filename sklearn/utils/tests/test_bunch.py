import warnings

import numpy as np
import pytest

from sklearn.utils import Bunch


def test_bunch_attribute_deprecation():
    """Check that bunch raises deprecation message with `__getattr__`."""
    bunch = Bunch()
    values = np.asarray([1, 2, 3])
    msg = (
        "Key: 'values', is deprecated in 1.3 and will be "
        "removed in 1.5. Please use 'grid_values' instead"
    )
    bunch._set_deprecated(
        values, new_key="grid_values", deprecated_key="values", warning_message=msg
    )

    with warnings.catch_warnings():
        # Does not warn for "grid_values"
        warnings.simplefilter("error")
        v = bunch["grid_values"]

    assert v is values

    with pytest.warns(FutureWarning, match=msg):
        # Warns for "values"
        v = bunch["values"]

    assert v is values


def test_bunch_copy():
    """Check that `copy` returns a `Bunch` and not a plain `dict`.

    Non-regression test for:
    https://github.com/scikit-learn/scikit-learn/issues/27930
    """
    values = np.asarray([1, 2, 3])
    bunch = Bunch(message="hello", data=values)
    bunch_copy = bunch.copy()

    assert isinstance(bunch_copy, Bunch)
    assert bunch_copy == bunch
    assert bunch_copy.message == "hello"
    # the copy is shallow
    assert bunch_copy.data is values

    bunch_copy.message = "goodbye"
    assert bunch.message == "hello"


def test_bunch_copy_keeps_deprecated_keys():
    """Check that deprecated keys still warn on a copied bunch."""
    bunch = Bunch()
    msg = "Key: 'values', is deprecated"
    bunch._set_deprecated(
        1, new_key="grid_values", deprecated_key="values", warning_message=msg
    )

    with warnings.catch_warnings():
        # Copying does not trigger the deprecation warning
        warnings.simplefilter("error")
        bunch_copy = bunch.copy()

    with pytest.warns(FutureWarning, match=msg):
        bunch_copy["values"]
