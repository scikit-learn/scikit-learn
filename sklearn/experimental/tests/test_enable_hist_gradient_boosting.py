"""Tests for making sure experimental imports work as expected."""

import textwrap

import pytest

from sklearn.utils._testing import assert_run_python_script_without_output
from sklearn.utils.fixes import _IS_WASM


# TODO(1.12): remove this test
@pytest.mark.xfail(_IS_WASM, reason="cannot start subprocess")
def test_import_raises_warning():
    code = """
    import pytest
    with pytest.warns(
        FutureWarning,
        match="deprecated since version 1.10 and will be removed in 1.12",
    ):
        from sklearn.experimental import enable_hist_gradient_boosting  # noqa
    """
    pattern = "deprecated since version 1.10 and will be removed in 1.12"
    assert_run_python_script_without_output(textwrap.dedent(code), pattern=pattern)
