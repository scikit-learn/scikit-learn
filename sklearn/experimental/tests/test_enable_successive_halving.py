"""Tests for making sure experimental imports work as expected."""

import textwrap

import pytest

from sklearn.utils._testing import assert_run_python_script_without_output
from sklearn.utils.fixes import _IS_WASM


# TODO(1.14): remove this test
@pytest.mark.xfail(_IS_WASM, reason="cannot start subprocess")
def test_import_raises_warning():
    code = """
    import pytest
    with pytest.warns(UserWarning, match="it is not needed to import"):
        from sklearn.experimental import enable_halving_search_cv  # noqa
    """
    pattern = "it is not needed to import"
    assert_run_python_script_without_output(textwrap.dedent(code), pattern=pattern)
