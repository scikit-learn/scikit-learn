import warnings

import pytest


def fff():
    warnings.warn("this is another warning", UserWarning)


def test_func():

    with pytest.warns(UserWarning, match=r"this is .* warning"):
        warnings.warn("this is a warning", UserWarning)
        fff()
