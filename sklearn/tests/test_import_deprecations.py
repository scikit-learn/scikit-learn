import textwrap

import pytest

from sklearn.utils._testing import assert_run_python_script
from sklearn._build_utils.deprecated_modules import _DEPRECATED_MODULES


# We are deprecating importing anything that isn't in an __init__ file and
# remaming most file.py into _file.py.
# This test makes sure imports are still possible but deprecated, with the
# appropriate error message.


@pytest.mark.parametrize('deprecated_path, importee', [
    (deprecated_path, importee)
    for _, deprecated_path, _, importee in _DEPRECATED_MODULES
])
def test_import_is_deprecated(deprecated_path, importee):
    # Make sure that "from deprecated_path import importee" is still possible
    # but raises a warning
    # We only need one entry per file, no need to check multiple imports from
    # the same file.

    # TODO: remove in 0.24

    # Special case for:
    # https://github.com/scikit-learn/scikit-learn/issues/15842
    if deprecated_path in ("sklearn.decomposition.dict_learning",
                           "sklearn.inspection.partial_dependence"):
        pytest.skip("No warning can be raised for " + deprecated_path)

    expected_message = (
        "The {deprecated_path} module is  deprecated in version "
        "0.22 and will be removed in version 0.24. "
        "The corresponding classes / functions "
        "should instead be imported from .*. "
        "Anything that cannot be imported from .* is now "
        "part of the private API."
    ).format(deprecated_path=deprecated_path)

    script = """
    import pytest

    with pytest.warns(FutureWarning,
                      match="{expected_message}"):
        from {deprecated_path} import {importee}
    """.format(
        expected_message=expected_message,
        deprecated_path=deprecated_path,
        importee=importee
    )
    assert_run_python_script(textwrap.dedent(script))
