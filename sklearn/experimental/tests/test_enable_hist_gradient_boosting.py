"""Tests for making sure experimental imports work as expected."""

import pytest


def test_import_raises_warning():
    with pytest.warns(UserWarning, match="It is not needed to import"):
        from sklearn.experimental import enable_hist_gradient_boosting
